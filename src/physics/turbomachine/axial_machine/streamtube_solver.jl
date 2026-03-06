using ....Utility: bracket_bisect_roots

function _mass_flow_invariant(
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
    f = nu_x -> mu - _mass_flow_invariant(gamma, pi, A, tau, nu_x, nu_theta)
    roots = bracket_bisect_roots(
        f,
        (x_lo, x_hi);
        n_scan=201,
        root_tol=1e-10,
        max_bisect_iters=80,
        dedupe_atol=1e-8,
    )
    isempty(roots) && return (converged=false, nu_x=NaN)
    nu_x = prefer == :high ? last(roots) : first(roots)
    return (converged=true, nu_x=nu_x)
end

function _invalid_streamtube_result(n_rows::Int; stall::Bool, choke::Bool, mu::Float64=NaN)
    n_stations = n_rows + 1
    return (
        PR=NaN,
        eta=NaN,
        stall=stall,
        choke=choke,
        valid=false,
        mu=mu,
        tau=fill(NaN, n_stations),
        pi=fill(NaN, n_stations),
        nu_theta=fill(NaN, n_stations),
        nu_x=fill(NaN, n_stations),
        stall_row=falses(n_rows),
        choke_row=falses(n_rows),
        valid_row=trues(n_rows),
    )
end

function _advance_row!(
    model::AxialMachineModel,
    row::AxialRow,
    aero::RotorAeroModel,
    row_radius::Float64,
    station_in::Int,
    station_out::Int,
    mu::Float64,
    tau::Vector{Float64},
    pi::Vector{Float64},
    nu_theta::Vector{Float64},
    nu_x::Vector{Float64},
    m_tip::Float64,
    prefer_root::Symbol,
)
    row.kind == :rotor || error("_advance_row! rotor method requires row.kind == :rotor")
    nu_u = row.speed_ratio_to_ref * m_tip * row_radius / model.r_tip_ref
    aero_out = row_aero(aero, nu_x[station_in], nu_theta[station_in], nu_u)
    stall = (aero_out.stall_margin <= 0) || !aero_out.valid
    aero_out.valid || return (converged=false, choke=false, stall=stall)

    gamma_ratio = model.gamma / (model.gamma - 1)
    mu_residual = function (nu_x_out)
        nu_theta_out = nu_u + nu_x_out * aero_out.k_theta_exit
        tau_out = tau[station_in] + (model.gamma - 1) * nu_u * (nu_theta_out - nu_theta[station_in])
        tau_out > 0 || return NaN
        pi_out = pi[station_in] * (tau_out / tau[station_in])^gamma_ratio *
                 exp(-gamma_ratio * aero_out.delta_s_hat)
        mu_at_nu_x = _mass_flow_invariant(
            model.gamma,
            pi_out,
            station_area(model, station_out),
            tau_out,
            nu_x_out,
            nu_theta_out,
        )
        return mu - mu_at_nu_x
    end

    roots = bracket_bisect_roots(
        mu_residual,
        (1e-10, 2.5);
        n_scan=201,
        root_tol=1e-10,
        max_bisect_iters=80,
        dedupe_atol=1e-8,
    )
    isempty(roots) && return (converged=false, choke=true, stall=stall)
    nu_x_out = prefer_root == :high ? last(roots) : first(roots)

    nu_theta_out = nu_u + nu_x_out * aero_out.k_theta_exit
    tau_out = tau[station_in] + (model.gamma - 1) * nu_u * (nu_theta_out - nu_theta[station_in])
    tau_out > 0 || return (converged=false, choke=false, stall=stall)

    pi_out = pi[station_in] *
             (tau_out / tau[station_in])^gamma_ratio *
             exp(-gamma_ratio * aero_out.delta_s_hat)

    nu_theta[station_out] = nu_theta_out
    tau[station_out] = tau_out
    pi[station_out] = pi_out
    nu_x[station_out] = nu_x_out
    return (converged=true, choke=false, stall=stall)
end

function _advance_row!(
    model::AxialMachineModel,
    row::AxialRow,
    aero::StatorAeroModel,
    _row_radius::Float64,
    station_in::Int,
    station_out::Int,
    mu::Float64,
    tau::Vector{Float64},
    pi::Vector{Float64},
    nu_theta::Vector{Float64},
    nu_x::Vector{Float64},
    _m_tip::Float64,
    prefer_root::Symbol,
)
    row.kind == :stator || error("_advance_row! stator method requires row.kind == :stator")
    aero_out = row_aero(aero, nu_x[station_in], nu_theta[station_in], 0.0)
    stall = (aero_out.stall_margin <= 0) || !aero_out.valid
    aero_out.valid || return (converged=false, choke=false, stall=stall)

    gamma_ratio = model.gamma / (model.gamma - 1)
    mu_residual = function (nu_x_out)
        nu_theta_out = nu_x_out * aero_out.k_theta_exit
        tau_out = tau[station_in]
        pi_out = pi[station_in] * exp(-gamma_ratio * aero_out.delta_s_hat)
        mu_at_nu_x = _mass_flow_invariant(
            model.gamma,
            pi_out,
            station_area(model, station_out),
            tau_out,
            nu_x_out,
            nu_theta_out,
        )
        return mu - mu_at_nu_x
    end

    roots = bracket_bisect_roots(
        mu_residual,
        (1e-10, 2.5);
        n_scan=201,
        root_tol=1e-10,
        max_bisect_iters=80,
        dedupe_atol=1e-8,
    )
    isempty(roots) && return (converged=false, choke=true, stall=stall)
    nu_x_out = prefer_root == :high ? last(roots) : first(roots)

    nu_theta[station_out] = nu_x_out * aero_out.k_theta_exit
    tau[station_out] = tau[station_in]
    pi[station_out] = pi[station_in] * exp(-gamma_ratio * aero_out.delta_s_hat)
    nu_x[station_out] = nu_x_out
    return (converged=true, choke=false, stall=stall)
end

function _nu_u_inlet_reference(
    model::AxialMachineModel,
    streamtube_radii::AbstractVector{<:Real},
    m_tip::Float64,
)
    idx_ref = model.first_rotor_index
    row_ref = model.rows[idx_ref]
    return row_ref.speed_ratio_to_ref * m_tip * Float64(streamtube_radii[idx_ref]) / model.r_tip_ref
end

"""
    streamtube_solve(model, streamtube_radii, m_tip, nu_x_inlet, nu_theta_inlet; prefer_root=:low)

Run the axial row-marching solve in non-dimensional coordinates.
"""
function streamtube_solve(
    model::AxialMachineModel,
    streamtube_radii::AbstractVector{<:Real},
    m_tip::Real,
    nu_x_inlet::Real,
    nu_theta_inlet::Real;
    prefer_root::Symbol=:low,
)
    m_tip_f = Float64(m_tip)
    nu_x_inlet_f = Float64(nu_x_inlet)
    n_rows = length(model.rows)
    n_stations = n_rows + 1
    length(streamtube_radii) == n_rows ||
        error("streamtube_radii length must match number of rows")
    radii = Float64.(streamtube_radii)
    for (k, row) in pairs(model.rows)
        row.r_hub <= radii[k] <= row.r_tip ||
            error("streamtube_radii[$k]=$(radii[k]) must lie in [r_hub, r_tip]=[$(row.r_hub), $(row.r_tip)]")
    end
    nu_x_inlet_f > 0 || error("nu_x_inlet must be > 0")

    tau = fill(NaN, n_stations)
    pi = fill(NaN, n_stations)
    nu_theta = fill(NaN, n_stations)
    nu_x = fill(NaN, n_stations)
    stall_row = falses(n_rows)
    choke_row = falses(n_rows)
    valid_row = trues(n_rows)

    tau[1] = 1.0
    pi[1] = 1.0
    nu_theta[1] = Float64(nu_theta_inlet)
    nu_x[1] = nu_x_inlet_f
    mu = _mass_flow_invariant(model.gamma, 1.0, station_area(model, 1), 1.0, nu_x[1], nu_theta[1])
    isfinite(mu) || return _invalid_streamtube_result(n_rows; stall=true, choke=true, mu=NaN)

    for k in 1:n_rows
        row = model.rows[k]
        station_in = k
        station_out = k + 1

        inlet = _solve_station_nu_x(
            model.gamma,
            mu,
            pi[station_in],
            station_area(model, station_in),
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

        row_step = _advance_row!(
            model,
            row,
            row.aero,
            radii[k],
            station_in,
            station_out,
            mu,
            tau,
            pi,
            nu_theta,
            nu_x,
            m_tip_f,
            prefer_root,
        )
        stall_row[k] = row_step.stall
        if !row_step.converged
            choke_row[k] = row_step.choke
            valid_row[k] = false
            break
        end
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

"""
    streamtube_solve_with_phi(model, m_tip, phi_in; streamtube_radii=meanline_radii(model), nu_theta_inlet=0.0, prefer_root=:low)

Phi-facing convenience wrapper around the nu_x-native solver.
"""
function streamtube_solve_with_phi(
    model::AxialMachineModel,
    m_tip::Real,
    phi_in::Real;
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    nu_theta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
)
    m_tip_f = Float64(m_tip)
    phi_in_f = Float64(phi_in)
    nu_u_ref = _nu_u_inlet_reference(model, streamtube_radii, m_tip_f)
    abs(nu_u_ref) > 0 || return _invalid_streamtube_result(length(model.rows); stall=true, choke=true)
    nu_x_inlet = phi_in_f * abs(nu_u_ref)
    return streamtube_solve(
        model,
        streamtube_radii,
        m_tip_f,
        nu_x_inlet,
        Float64(nu_theta_inlet);
        prefer_root=prefer_root,
    )
end

"""
    sample_streamtube_solve(model, speed_grid, flow_grid; flow_min=nothing, flow_max=nothing, streamtube_radii=meanline_radii(model), nu_theta_inlet=0.0, prefer_root=:low, is_feasible)

Sample the phi-facing streamtube solver over a Cartesian grid of speed/flow coordinates.
"""
function sample_streamtube_solve(
    model::AxialMachineModel,
    speed_grid::AbstractVector{<:Real},
    flow_grid::AbstractVector{<:Real};
    flow_min::Union{Nothing,AbstractVector{<:Real}}=nothing,
    flow_max::Union{Nothing,AbstractVector{<:Real}}=nothing,
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    nu_theta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
    is_feasible::Function=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta)),
)
    length(speed_grid) >= 1 || error("speed_grid must be non-empty")
    length(flow_grid) >= 1 || error("flow_grid must be non-empty")

    has_limits = !isnothing(flow_min) || !isnothing(flow_max)
    if has_limits
        isnothing(flow_min) && error("flow_min must be provided when using flow limits")
        isnothing(flow_max) && error("flow_max must be provided when using flow limits")
        length(flow_min) == length(speed_grid) || error("flow_min length must match speed_grid")
        length(flow_max) == length(speed_grid) || error("flow_max length must match speed_grid")
    end

    pr_table = Matrix{Float64}(undef, length(speed_grid), length(flow_grid))
    eta_table = Matrix{Float64}(undef, length(speed_grid), length(flow_grid))

    for (i, speed_raw) in pairs(speed_grid)
        speed = Float64(speed_raw)
        for (j, flow_raw) in pairs(flow_grid)
            flow = Float64(flow_raw)
            if has_limits
                flow = clamp(flow, Float64(flow_min[i]), Float64(flow_max[i]))
            end
            vals = streamtube_solve_with_phi(
                model,
                speed,
                flow;
                streamtube_radii=streamtube_radii,
                nu_theta_inlet=nu_theta_inlet,
                prefer_root=prefer_root,
            )
            is_feasible(vals) ||
                error("streamtube sampling produced invalid value at speed=$(speed), flow=$(flow)")
            pr_table[i, j] = vals.PR
            eta_table[i, j] = vals.eta
        end
    end

    return (pr_table=pr_table, eta_table=eta_table)
end
