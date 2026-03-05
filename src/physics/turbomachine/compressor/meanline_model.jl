"""
Compressor meanline row-marching model.
"""

using TOML
import ....Utility: write_toml, read_toml

"""
Abstract aerodynamic row model used by compressor meanline marching.
"""
abstract type AbstractRowAeroModel end

"""
Simple rotor aerodynamic closure.

Parameter definitions and tuning notes:
`docs/compressor_meanline_row_aero_parameters.md`
"""
Base.@kwdef struct RotorAeroModel{T<:Real} <: AbstractRowAeroModel
    beta_ref::T = T(-0.55)                  # reference relative inlet angle [rad]
    beta_incidence_sensitivity::T = T(0.75) # exit-angle response to incidence [-]
    loss_base::T = T(0.010)                 # baseline entropy rise, Delta s0 / cp [-]
    loss_incidence::T = T(0.18)             # quadratic incidence-loss gain [-]
    stall_incidence_limit::T = T(0.32)      # |incidence| stall threshold [rad]
    k_theta_min::T = T(-2.5)                # tan(beta_out) lower clamp [-]
    k_theta_max::T = T(1.5)                 # tan(beta_out) upper clamp [-]
end

"""
Simple stator aerodynamic closure.

Parameter definitions and tuning notes:
`docs/compressor_meanline_row_aero_parameters.md`
"""
Base.@kwdef struct StatorAeroModel{T<:Real} <: AbstractRowAeroModel
    alpha_ref::T = T(0.45)                   # reference absolute inlet angle [rad]
    alpha_incidence_sensitivity::T = T(0.85) # exit-angle response to incidence [-]
    loss_base::T = T(0.006)                  # baseline entropy rise, Delta s0 / cp [-]
    loss_incidence::T = T(0.12)              # quadratic incidence-loss gain [-]
    stall_incidence_limit::T = T(0.30)       # |incidence| stall threshold [rad]
    k_theta_min::T = T(-1.2)                 # tan(alpha_out) lower clamp [-]
    k_theta_max::T = T(2.5)                  # tan(alpha_out) upper clamp [-]
end

"""
Row aerodynamic input state.
"""
struct RowAeroInput{T}
    nu_x_in::T
    nu_theta_in::T
    nu_u::T
    tau_in::T
    pi_in::T
    mu::T
end

"""
Row aerodynamic output contract for the meanline marcher.
"""
struct RowAeroOutput{T}
    k_theta_exit::T
    delta_s_hat::T
    stall_margin::T
    valid::Bool
    diagnostics::NamedTuple
end

"""
Compressor row descriptor.
"""
struct CompressorRow{M<:AbstractRowAeroModel}
    kind::Symbol            # :rotor or :stator
    aero::M
    r_mean::Float64
    r_tip::Float64
    omega_sign::Int8        # +1/-1 for rotor, 0 for stator
end

function CompressorRow(
    kind::Symbol,
    aero::AbstractRowAeroModel,
    r_mean::Real,
    r_tip::Real,
    omega_sign::Integer,
)
    kind in (:rotor, :stator) || error("row kind must be :rotor or :stator")
    r_mean > 0 || error("row r_mean must be > 0")
    r_tip > 0 || error("row r_tip must be > 0")
    sign_i8 = Int8(omega_sign)
    if kind == :rotor
        sign_i8 in Int8[-1, 1] || error("rotor omega_sign must be -1 or +1")
        aero isa RotorAeroModel || error("rotor rows require RotorAeroModel")
    else
        sign_i8 == 0 || error("stator omega_sign must be 0")
        aero isa StatorAeroModel || error("stator rows require StatorAeroModel")
    end
    return CompressorRow{typeof(aero)}(kind, aero, Float64(r_mean), Float64(r_tip), sign_i8)
end

"""
Compressor meanline model with row sequence and area schedule.
"""
struct CompressorMeanlineModel <: AbstractCompressorPerformanceMap
    gamma::Float64
    gas_constant::Float64
    A_ref::Float64
    A_station::Vector{Float64}
    rows::Vector{CompressorRow}
    nu_theta_station1::Float64
    m_tip_bounds::Tuple{Float64,Float64}
    phi_in_bounds::Tuple{Float64,Float64}
end

function CompressorMeanlineModel(
    gamma::Real,
    gas_constant::Real,
    A_ref::Real,
    A_station::Vector{<:Real},
    rows::Vector{CompressorRow},
    nu_theta_station1::Real,
    m_tip_bounds::Tuple{<:Real,<:Real},
    phi_in_bounds::Tuple{<:Real,<:Real},
)
    gamma > 1 || error("gamma must be > 1")
    gas_constant > 0 || error("gas_constant must be > 0")
    A_ref > 0 || error("A_ref must be > 0")
    isempty(rows) && error("rows must not be empty")
    length(A_station) == length(rows) + 1 || error("A_station length must be n_rows + 1")
    all(a -> a > 0, A_station) || error("A_station entries must be > 0")
    m_lo, m_hi = Float64(m_tip_bounds[1]), Float64(m_tip_bounds[2])
    phi_lo, phi_hi = Float64(phi_in_bounds[1]), Float64(phi_in_bounds[2])
    m_hi > m_lo > 0 || error("m_tip_bounds must satisfy 0 < lo < hi")
    phi_hi > phi_lo > 0 || error("phi_in_bounds must satisfy 0 < lo < hi")
    return CompressorMeanlineModel(
        Float64(gamma),
        Float64(gas_constant),
        Float64(A_ref),
        Float64.(A_station),
        rows,
        Float64(nu_theta_station1),
        (m_lo, m_hi),
        (phi_lo, phi_hi),
    )
end

function _first_rotor_index(model::CompressorMeanlineModel)
    idx = findfirst(row -> row.kind == :rotor, model.rows)
    return isnothing(idx) ? 1 : idx
end

function _station_rhs(
    gamma::Real,
    pi::Real,
    A::Real,
    tau::Real,
    nu_x::Real,
    nu_theta::Real,
)
    tau > 0 || return NaN
    term = 1 - ((gamma - 1) / (2 * tau)) * (nu_x^2 + nu_theta^2)
    term > 0 || return NaN
    return pi * A * (nu_x / tau) * term^(1 / (gamma - 1))
end

function _bisect_root(
    f::Function,
    a::Float64,
    b::Float64;
    tol::Float64=1e-10,
    max_iters::Int=80,
)
    fa = f(a)
    fb = f(b)
    isfinite(fa) && isfinite(fb) || return NaN
    fa * fb <= 0 || return NaN
    lo = a
    hi = b
    flo = fa
    mid = 0.5 * (lo + hi)
    for _ in 1:max_iters
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        isfinite(fmid) || return NaN
        if abs(fmid) <= tol || abs(hi - lo) <= tol
            return mid
        end
        if flo * fmid <= 0
            hi = mid
        else
            lo = mid
            flo = fmid
        end
    end
    return mid
end

function _find_scalar_roots(
    f::Function,
    x_lo::Float64,
    x_hi::Float64;
    n_scan::Int=201,
    tol::Float64=1e-10,
)
    x_hi > x_lo || return Float64[]
    grid = collect(range(x_lo, x_hi, length=n_scan))
    vals = [f(x) for x in grid]
    roots = Float64[]
    for i in eachindex(grid)
        fi = vals[i]
        if isfinite(fi) && abs(fi) <= tol
            push!(roots, grid[i])
        end
    end
    for i in 1:(length(grid) - 1)
        f1 = vals[i]
        f2 = vals[i + 1]
        if !(isfinite(f1) && isfinite(f2))
            continue
        end
        if f1 * f2 < 0
            r = _bisect_root(f, grid[i], grid[i + 1]; tol=tol)
            isfinite(r) && push!(roots, r)
        end
    end
    isempty(roots) && return roots
    sort!(roots)
    out = Float64[roots[1]]
    for i in 2:length(roots)
        abs(roots[i] - out[end]) > 1e-8 && push!(out, roots[i])
    end
    return out
end

function _solve_station_nu_x(
    gamma::Float64,
    mu::Float64,
    pi::Float64,
    A::Float64,
    tau::Float64,
    nu_theta::Float64;
    prefer::Symbol=:low,
)
    tau > 0 || return (converged=false, nu_x=NaN)
    term = (2 * tau / (gamma - 1)) - nu_theta^2
    term > 0 || return (converged=false, nu_x=NaN)
    x_hi = sqrt(term) * (1 - 1e-8)
    x_lo = 1e-10
    f = nu_x -> mu - _station_rhs(gamma, pi, A, tau, nu_x, nu_theta)
    roots = _find_scalar_roots(f, x_lo, x_hi)
    isempty(roots) && return (converged=false, nu_x=NaN)
    nu_x = prefer == :high ? last(roots) : first(roots)
    return (converged=true, nu_x=nu_x)
end

function row_aero(
    model::RotorAeroModel{T},
    row::CompressorRow,
    input::RowAeroInput{U},
) where {T<:Real,U<:Real}
    beta_in = atan(input.nu_theta_in - input.nu_u, input.nu_x_in)
    incidence = beta_in - model.beta_ref
    beta_out = model.beta_ref - model.beta_incidence_sensitivity * incidence
    k_theta = clamp(tan(beta_out), model.k_theta_min, model.k_theta_max)
    delta_s_hat = model.loss_base + model.loss_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(k_theta) && isfinite(delta_s_hat)
    return RowAeroOutput(
        k_theta,
        max(delta_s_hat, zero(delta_s_hat)),
        stall_margin,
        valid,
        (incidence=incidence, beta_in=beta_in, beta_out=beta_out),
    )
end

function row_aero(
    model::StatorAeroModel{T},
    row::CompressorRow,
    input::RowAeroInput{U},
) where {T<:Real,U<:Real}
    alpha_in = atan(input.nu_theta_in, input.nu_x_in)
    incidence = alpha_in - model.alpha_ref
    alpha_out = model.alpha_ref - model.alpha_incidence_sensitivity * incidence
    k_theta = clamp(tan(alpha_out), model.k_theta_min, model.k_theta_max)
    delta_s_hat = model.loss_base + model.loss_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(k_theta) && isfinite(delta_s_hat)
    return RowAeroOutput(
        k_theta,
        max(delta_s_hat, zero(delta_s_hat)),
        stall_margin,
        valid,
        (incidence=incidence, alpha_in=alpha_in, alpha_out=alpha_out),
    )
end

function _nu_u_for_row(
    model::CompressorMeanlineModel,
    row::CompressorRow,
    m_tip::Float64,
)
    idx_ref = _first_rotor_index(model)
    r_tip_ref = model.rows[idx_ref].r_tip
    if row.kind == :rotor
        return Float64(row.omega_sign) * m_tip * row.r_mean / r_tip_ref
    else
        return 0.0
    end
end

function _meanline_performance_map(
    model::CompressorMeanlineModel,
    m_tip::Real,
    phi_in::Real;
    prefer_root::Symbol=:low,
)
    m_tip_f = Float64(m_tip)
    phi_in_f = Float64(phi_in)
    n_rows = length(model.rows)
    n_stations = n_rows + 1

    tau = fill(NaN, n_stations)
    pi = fill(NaN, n_stations)
    nu_theta = fill(NaN, n_stations)
    nu_x = fill(NaN, n_stations)
    stall_row = falses(n_rows)
    choke_row = falses(n_rows)
    valid_row = trues(n_rows)

    tau[1] = 1.0
    pi[1] = 1.0
    nu_theta[1] = model.nu_theta_station1
    nu_u1 = _nu_u_for_row(model, model.rows[_first_rotor_index(model)], m_tip_f)
    abs(nu_u1) > 0 || return (
        PR=NaN,
        eta=NaN,
        stall=true,
        choke=true,
        valid=false,
        mu=NaN,
        tau=tau,
        pi=pi,
        nu_theta=nu_theta,
        nu_x=nu_x,
        stall_row=stall_row,
        choke_row=choke_row,
        valid_row=valid_row,
    )
    nu_x[1] = phi_in_f * abs(nu_u1)
    mu = _station_rhs(model.gamma, 1.0, model.A_station[1], 1.0, nu_x[1], nu_theta[1])
    isfinite(mu) || return (
        PR=NaN,
        eta=NaN,
        stall=true,
        choke=true,
        valid=false,
        mu=NaN,
        tau=tau,
        pi=pi,
        nu_theta=nu_theta,
        nu_x=nu_x,
        stall_row=stall_row,
        choke_row=choke_row,
        valid_row=valid_row,
    )

    for k in 1:n_rows
        row = model.rows[k]
        station_in = k
        station_out = k + 1

        inlet = _solve_station_nu_x(
            model.gamma,
            mu,
            pi[station_in],
            model.A_station[station_in],
            tau[station_in],
            nu_theta[station_in];
            prefer=prefer_root,
        )
        if !inlet.converged
            choke_row[k] = true
            valid_row[k] = false
            break
        end
        nu_x[station_in] = inlet.nu_x

        nu_u = _nu_u_for_row(model, row, m_tip_f)
        aero_in = RowAeroInput(
            nu_x[station_in],
            nu_theta[station_in],
            nu_u,
            tau[station_in],
            pi[station_in],
            mu,
        )
        aero_out = row_aero(row.aero, row, aero_in)
        stall_row[k] = (aero_out.stall_margin <= 0) || !aero_out.valid
        if !aero_out.valid
            valid_row[k] = false
            break
        end

        gamma_ratio = model.gamma / (model.gamma - 1)
        f_exit = function (nu_x_out)
            if row.kind == :rotor
                nu_theta_out = nu_u + nu_x_out * aero_out.k_theta_exit
                tau_out = tau[station_in] + (model.gamma - 1) * nu_u * (nu_theta_out - nu_theta[station_in])
                tau_out > 0 || return NaN
                pi_out = pi[station_in] * (tau_out / tau[station_in])^gamma_ratio *
                         exp(-gamma_ratio * aero_out.delta_s_hat)
                rhs = _station_rhs(
                    model.gamma,
                    pi_out,
                    model.A_station[station_out],
                    tau_out,
                    nu_x_out,
                    nu_theta_out,
                )
                return mu - rhs
            else
                nu_theta_out = nu_x_out * aero_out.k_theta_exit
                tau_out = tau[station_in]
                pi_out = pi[station_in] * exp(-gamma_ratio * aero_out.delta_s_hat)
                rhs = _station_rhs(
                    model.gamma,
                    pi_out,
                    model.A_station[station_out],
                    tau_out,
                    nu_x_out,
                    nu_theta_out,
                )
                return mu - rhs
            end
        end

        # Broad positive nu_x search range; physical admissibility is enforced
        # inside residual evaluation via _station_rhs static-temperature checks.
        roots = _find_scalar_roots(f_exit, 1e-10, 2.5)
        if isempty(roots)
            choke_row[k] = true
            valid_row[k] = false
            break
        end
        nu_x_out = prefer_root == :high ? last(roots) : first(roots)

        if row.kind == :rotor
            nu_theta[station_out] = nu_u + nu_x_out * aero_out.k_theta_exit
            tau[station_out] = tau[station_in] + (model.gamma - 1) * nu_u *
                               (nu_theta[station_out] - nu_theta[station_in])
            tau[station_out] > 0 || begin
                valid_row[k] = false
                break
            end
            pi[station_out] = pi[station_in] *
                              (tau[station_out] / tau[station_in])^gamma_ratio *
                              exp(-gamma_ratio * aero_out.delta_s_hat)
        else
            nu_theta[station_out] = nu_x_out * aero_out.k_theta_exit
            tau[station_out] = tau[station_in]
            pi[station_out] = pi[station_in] * exp(-gamma_ratio * aero_out.delta_s_hat)
        end
        nu_x[station_out] = nu_x_out
    end

    model_valid = all(valid_row)
    model_choke = any(choke_row)
    model_stall = any(stall_row)
    if !model_valid
        return (
            PR=NaN,
            eta=NaN,
            stall=model_stall,
            choke=model_choke,
            valid=false,
            mu=mu,
            tau=tau,
            pi=pi,
            nu_theta=nu_theta,
            nu_x=nu_x,
            stall_row=stall_row,
            choke_row=choke_row,
            valid_row=valid_row,
        )
    end

    tau_out = tau[end]
    pi_out = pi[end]
    # Allow both compressor mode (PR > 1) and windmilling mode (PR <= 1).
    # eta_tt is evaluated whenever mathematically well-defined.
    eta_tt = if (tau_out > 0) && (pi_out > 0) && (abs(tau_out - 1) > 1e-12)
        eta = (pi_out^((model.gamma - 1) / model.gamma) - 1) / (tau_out - 1)
        isfinite(eta) ? eta : NaN
    else
        NaN
    end
    return (
        PR=pi_out,
        eta=eta_tt,
        stall=model_stall,
        choke=model_choke,
        valid=model_valid && isfinite(pi_out) && (pi_out > 0) && isfinite(eta_tt),
        mu=mu,
        tau=tau,
        pi=pi,
        nu_theta=nu_theta,
        nu_x=nu_x,
        stall_row=stall_row,
        choke_row=choke_row,
        valid_row=valid_row,
    )
end

function compressor_performance_map_from_stagnation(
    model::CompressorMeanlineModel,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    Tt_in > 0 || error("Tt_in must be > 0")
    Pt_in > 0 || error("Pt_in must be > 0")
    rho0 = Pt_in / (model.gas_constant * Tt_in)
    a0_in = sqrt(model.gamma * model.gas_constant * Tt_in)
    idx_ref = _first_rotor_index(model)
    r_tip_1 = model.rows[idx_ref].r_tip
    r_mean_1 = model.rows[idx_ref].r_mean
    A_phys_1 = model.A_ref * model.A_station[1]
    m_tip = omega * r_tip_1 / a0_in
    U_m1 = abs(omega) * r_mean_1
    U_m1 > 0 || error("omega must be non-zero for meanline compressor map evaluation")
    Vx_1 = mdot / (rho0 * A_phys_1)
    phi_in = Vx_1 / U_m1
    vals = _meanline_performance_map(model, m_tip, phi_in)
    return (
        PR=vals.PR,
        eta=vals.eta,
        speed_coord=m_tip,
        flow_coord=phi_in,
        stall=vals.stall,
        choke=vals.choke,
        valid=vals.valid,
    )
end

function performance_map_domain(
    model::CompressorMeanlineModel,
    Tt_in::Real,
    Pt_in::Real,
)
    Tt_in > 0 || error("Tt_in must be > 0")
    Pt_in > 0 || error("Pt_in must be > 0")
    idx_ref = _first_rotor_index(model)
    r_tip_1 = model.rows[idx_ref].r_tip
    r_mean_1 = model.rows[idx_ref].r_mean
    A_phys_1 = model.A_ref * model.A_station[1]
    rho0 = Pt_in / (model.gas_constant * Tt_in)
    a0 = sqrt(model.gamma * model.gas_constant * Tt_in)
    m_lo, m_hi = model.m_tip_bounds
    phi_lo, phi_hi = model.phi_in_bounds
    omega_lo = m_lo * a0 / r_tip_1
    omega_hi = m_hi * a0 / r_tip_1

    mdot_from = (omega, phi) -> begin
        U_m1 = abs(omega) * r_mean_1
        return phi * U_m1 * rho0 * A_phys_1
    end

    mdot_vals = Float64[
        mdot_from(omega_lo, phi_lo),
        mdot_from(omega_lo, phi_hi),
        mdot_from(omega_hi, phi_lo),
        mdot_from(omega_hi, phi_hi),
    ]
    return (
        omega=(omega_lo, omega_hi),
        mdot=(minimum(mdot_vals), maximum(mdot_vals)),
        mdot_flow_range=(
            surge=(omega -> mdot_from(omega, phi_lo)),
            choke=(omega -> mdot_from(omega, phi_hi)),
        ),
    )
end

function _meanline_coords_to_operating_point(
    model::CompressorMeanlineModel,
    m_tip::Real,
    phi_in::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    idx_ref = _first_rotor_index(model)
    r_tip_1 = model.rows[idx_ref].r_tip
    r_mean_1 = model.rows[idx_ref].r_mean
    A_phys_1 = model.A_ref * model.A_station[1]
    a0_in = sqrt(model.gamma * model.gas_constant * Tt_in)
    rho0 = Pt_in / (model.gas_constant * Tt_in)
    omega = m_tip * a0_in / r_tip_1
    U_m1 = abs(omega) * r_mean_1
    mdot = phi_in * U_m1 * rho0 * A_phys_1
    return (omega=omega, mdot=mdot)
end

"""
Tabulate a meanline compressor model onto a non-dimensional performance-map grid.
"""
function tabulate_compressor_meanline_model(
    model::CompressorMeanlineModel;
    Tt_in_ref::Real=288.15,
    Pt_in_ref::Real=101_325.0,
    n_speed::Int=31,
    n_flow::Int=41,
    m_tip_grid::Union{Nothing,Vector{<:Real}}=nothing,
    phi_in_grid::Union{Nothing,Vector{<:Real}}=nothing,
    interpolation::Symbol=:bilinear,
    boundary_resolution::Int=401,
)
    Tt_in_ref > 0 || error("Tt_in_ref must be > 0")
    Pt_in_ref > 0 || error("Pt_in_ref must be > 0")
    n_speed >= 2 || error("n_speed must be >= 2")
    n_flow >= 2 || error("n_flow must be >= 2")
    boundary_resolution >= 21 || error("boundary_resolution must be >= 21")
    interpolation in (:bilinear, :bicubic) ||
        error("interpolation must be :bilinear or :bicubic")

    m_lo, m_hi = model.m_tip_bounds
    phi_lo, phi_hi = model.phi_in_bounds
    m_grid = isnothing(m_tip_grid) ? collect(range(m_lo, m_hi, length=n_speed)) : Float64.(m_tip_grid)
    phi_grid = isnothing(phi_in_grid) ? collect(range(phi_lo, phi_hi, length=n_flow)) : Float64.(phi_in_grid)
    length(m_grid) >= 2 || error("m_tip_grid must have at least 2 points")
    length(phi_grid) >= 2 || error("phi_in_grid must have at least 2 points")
    issorted(m_grid) || error("m_tip_grid must be sorted ascending")
    issorted(phi_grid) || error("phi_in_grid must be sorted ascending")

    function eval_at(m_tip::Float64, phi::Float64)
        op = _meanline_coords_to_operating_point(model, m_tip, phi, Tt_in_ref, Pt_in_ref)
        return compressor_performance_map_from_stagnation(
            model,
            op.omega,
            op.mdot,
            Tt_in_ref,
            Pt_in_ref,
        )
    end

    phi_probe = collect(range(phi_lo, phi_hi, length=boundary_resolution))
    valid_speed_idx = Int[]
    phi_surge = Float64[]
    phi_choke = Float64[]
    for (i, m_tip) in pairs(m_grid)
        valid_phis = Float64[]
        for phi in phi_probe
            vals = eval_at(m_tip, phi)
            if vals.valid && isfinite(vals.PR) && isfinite(vals.eta)
                push!(valid_phis, phi)
            end
        end
        isempty(valid_phis) && continue
        push!(valid_speed_idx, i)
        push!(phi_surge, first(valid_phis))
        push!(phi_choke, last(valid_phis))
    end
    length(valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid speed lines in requested tabulation range")
    m_grid_valid = m_grid[valid_speed_idx]

    pr_table = Matrix{Float64}(undef, length(m_grid_valid), length(phi_grid))
    eta_table = Matrix{Float64}(undef, length(m_grid_valid), length(phi_grid))
    for (i, m_tip) in pairs(m_grid_valid)
        for (j, phi_raw) in pairs(phi_grid)
            phi = clamp(phi_raw, phi_surge[i], phi_choke[i])
            vals = eval_at(m_tip, phi)
            (vals.valid && isfinite(vals.PR) && isfinite(vals.eta)) ||
                error("meanline tabulation produced non-finite value at m_tip=$(m_tip), phi=$(phi)")
            pr_table[i, j] = vals.PR
            eta_table[i, j] = vals.eta
        end
    end

    idx_ref = _first_rotor_index(model)
    tip_radius_inlet = model.rows[idx_ref].r_tip
    mean_radius_inlet = model.rows[idx_ref].r_mean
    inlet_area = model.A_ref * model.A_station[1]
    return NonDimensionalTabulatedCompressorPerformanceMap(
        model.gamma,
        model.gas_constant,
        tip_radius_inlet,
        mean_radius_inlet,
        inlet_area,
        m_grid_valid,
        phi_grid,
        pr_table,
        eta_table;
        interpolation=interpolation,
        phi_surge=phi_surge,
        phi_choke=phi_choke,
    )
end

function _aero_to_toml_dict(model::RotorAeroModel)
    return Dict{String,Any}(
        "format" => "rotor_aero_model",
        "beta_ref" => Float64(model.beta_ref),
        "beta_incidence_sensitivity" => Float64(model.beta_incidence_sensitivity),
        "loss_base" => Float64(model.loss_base),
        "loss_incidence" => Float64(model.loss_incidence),
        "stall_incidence_limit" => Float64(model.stall_incidence_limit),
        "k_theta_min" => Float64(model.k_theta_min),
        "k_theta_max" => Float64(model.k_theta_max),
    )
end

function _aero_to_toml_dict(model::StatorAeroModel)
    return Dict{String,Any}(
        "format" => "stator_aero_model",
        "alpha_ref" => Float64(model.alpha_ref),
        "alpha_incidence_sensitivity" => Float64(model.alpha_incidence_sensitivity),
        "loss_base" => Float64(model.loss_base),
        "loss_incidence" => Float64(model.loss_incidence),
        "stall_incidence_limit" => Float64(model.stall_incidence_limit),
        "k_theta_min" => Float64(model.k_theta_min),
        "k_theta_max" => Float64(model.k_theta_max),
    )
end

function _aero_from_toml_dict(data::Dict{String,Any})
    haskey(data, "format") || error("aero model missing format")
    fmt = String(data["format"])
    if fmt == "rotor_aero_model"
        return RotorAeroModel{Float64}(
            beta_ref=Float64(data["beta_ref"]),
            beta_incidence_sensitivity=Float64(data["beta_incidence_sensitivity"]),
            loss_base=Float64(data["loss_base"]),
            loss_incidence=Float64(data["loss_incidence"]),
            stall_incidence_limit=Float64(data["stall_incidence_limit"]),
            k_theta_min=Float64(data["k_theta_min"]),
            k_theta_max=Float64(data["k_theta_max"]),
        )
    elseif fmt == "stator_aero_model"
        return StatorAeroModel{Float64}(
            alpha_ref=Float64(data["alpha_ref"]),
            alpha_incidence_sensitivity=Float64(data["alpha_incidence_sensitivity"]),
            loss_base=Float64(data["loss_base"]),
            loss_incidence=Float64(data["loss_incidence"]),
            stall_incidence_limit=Float64(data["stall_incidence_limit"]),
            k_theta_min=Float64(data["k_theta_min"]),
            k_theta_max=Float64(data["k_theta_max"]),
        )
    end
    error("unsupported aero model format=$(fmt)")
end

function _row_to_toml_dict(row::CompressorRow)
    return Dict{String,Any}(
        "kind" => String(row.kind),
        "r_mean" => row.r_mean,
        "r_tip" => row.r_tip,
        "omega_sign" => Int(row.omega_sign),
        "aero" => _aero_to_toml_dict(row.aero),
    )
end

function _row_from_toml_dict(data::Dict{String,Any})
    haskey(data, "kind") || error("row missing kind")
    haskey(data, "r_mean") || error("row missing r_mean")
    haskey(data, "r_tip") || error("row missing r_tip")
    haskey(data, "omega_sign") || error("row missing omega_sign")
    haskey(data, "aero") || error("row missing aero")
    kind = Symbol(String(data["kind"]))
    aero = _aero_from_toml_dict(data["aero"])
    return CompressorRow(
        kind,
        aero,
        Float64(data["r_mean"]),
        Float64(data["r_tip"]),
        Int(data["omega_sign"]),
    )
end

function write_toml(
    model::CompressorMeanlineModel,
    path::AbstractString;
    group::AbstractString="compressor_meanline_model",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "compressor_meanline_model"
    node["format_version"] = 1
    node["gamma"] = model.gamma
    node["gas_constant"] = model.gas_constant
    node["A_ref"] = model.A_ref
    node["A_station"] = model.A_station
    node["nu_theta_station1"] = model.nu_theta_station1
    node["m_tip_bounds"] = [model.m_tip_bounds[1], model.m_tip_bounds[2]]
    node["phi_in_bounds"] = [model.phi_in_bounds[1], model.phi_in_bounds[2]]
    node["rows"] = [_row_to_toml_dict(row) for row in model.rows]
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{CompressorMeanlineModel},
    path::AbstractString;
    group::AbstractString="compressor_meanline_model",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    for key in ("gamma", "gas_constant", "A_ref", "A_station", "rows")
        haskey(node, key) || error("missing TOML key $(key)")
    end
    nu_theta_station1 = haskey(node, "nu_theta_station1") ? Float64(node["nu_theta_station1"]) : 0.0
    m_tip_bounds = haskey(node, "m_tip_bounds") ? (Float64(node["m_tip_bounds"][1]), Float64(node["m_tip_bounds"][2])) : (0.5, 1.2)
    phi_in_bounds = haskey(node, "phi_in_bounds") ? (Float64(node["phi_in_bounds"][1]), Float64(node["phi_in_bounds"][2])) : (0.2, 0.9)
    rows = CompressorRow[_row_from_toml_dict(row_data) for row_data in node["rows"]]
    return CompressorMeanlineModel(
        Float64(node["gamma"]),
        Float64(node["gas_constant"]),
        Float64(node["A_ref"]),
        Float64.(node["A_station"]),
        rows,
        nu_theta_station1,
        m_tip_bounds,
        phi_in_bounds,
    )
end

"""
Demo compressor meanline model for development/testing.
"""
function demo_compressor_meanline_model()
    rows = CompressorRow[
        CompressorRow(
            :rotor,
            RotorAeroModel{Float64}(
                beta_ref=-0.55,
                beta_incidence_sensitivity=0.62,
                loss_base=0.0025,
                loss_incidence=0.045,
                stall_incidence_limit=0.36,
                k_theta_min=-2.0,
                k_theta_max=1.1,
            ),
            0.180,
            0.220,
            +1,
        ),
        CompressorRow(
            :stator,
            StatorAeroModel{Float64}(
                alpha_ref=0.45,
                alpha_incidence_sensitivity=0.70,
                loss_base=0.0018,
                loss_incidence=0.030,
                stall_incidence_limit=0.34,
                k_theta_min=-1.0,
                k_theta_max=2.0,
            ),
            0.180,
            0.220,
            0,
        ),
        CompressorRow(
            :rotor,
            RotorAeroModel{Float64}(
                beta_ref=-0.50,
                beta_incidence_sensitivity=0.60,
                loss_base=0.0030,
                loss_incidence=0.050,
                stall_incidence_limit=0.34,
                k_theta_min=-2.1,
                k_theta_max=1.1,
            ),
            0.180,
            0.220,
            +1,
        ),
        CompressorRow(
            :stator,
            StatorAeroModel{Float64}(
                alpha_ref=0.45,
                alpha_incidence_sensitivity=0.70,
                loss_base=0.0020,
                loss_incidence=0.032,
                stall_incidence_limit=0.34,
                k_theta_min=-1.0,
                k_theta_max=2.0,
            ),
            0.180,
            0.220,
            0,
        ),
    ]
    return CompressorMeanlineModel(
        1.4,
        287.05,
        0.060,
        [1.00, 0.98, 0.96, 0.95, 0.95],
        rows,
        0.0,
        (0.01, 1.10),
        (0.01, 0.95),
    )
end
