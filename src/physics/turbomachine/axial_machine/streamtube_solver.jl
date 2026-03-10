using ....Utility: bracket_bisect_roots

function _station_radius(model::AxialMachineModel, streamtube_radii::AbstractVector{<:Real}, station_index::Int)
    n_rows = length(model.rows)
    1 <= station_index <= (n_rows + 1) || error("station_index out of bounds")
    if station_index == 1
        return Float64(streamtube_radii[1])
    elseif station_index == n_rows + 1
        return Float64(streamtube_radii[end])
    end
    return 0.5 * (Float64(streamtube_radii[station_index - 1]) + Float64(streamtube_radii[station_index]))
end

function _build_station_data(
    model::AxialMachineModel,
    streamtube_radii::AbstractVector{<:Real},
    tau::Vector{Float64},
    pi::Vector{Float64},
    nu_theta::Vector{Float64},
    nu_x::Vector{Float64},
)
    n_stations = length(tau)
    return [
        (
            station_index=k,
            radius=_station_radius(model, streamtube_radii, k),
            area=station_area(model, k),
            tau=tau[k],
            pi=pi[k],
            nu_theta=nu_theta[k],
            nu_x=nu_x[k],
        ) for k in 1:n_stations
    ]
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
        station_data=NamedTuple[],
        row_data=NamedTuple[],
    )
end

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

function _solve_row_nu_x(
    model::AxialMachineModel,
    station_in::Int,
    station_out::Int,
    mu::Float64,
    tau::Vector{Float64},
    pi::Vector{Float64},
    nu_theta::Vector{Float64},
    nu_u::Float64,
    k_theta_exit::Float64,
    delta_s_hat::Float64,
    prefer_root::Symbol,
)
    gamma_ratio = model.gamma / (model.gamma - 1)
    mu_residual = function (nu_x_out)
        nu_theta_out = nu_u + nu_x_out * k_theta_exit
        tau_out = tau[station_in] + (model.gamma - 1) * nu_u * (nu_theta_out - nu_theta[station_in])
        tau_out > 0 || return NaN
        pi_out = pi[station_in] * (tau_out / tau[station_in])^gamma_ratio *
                 exp(-gamma_ratio * delta_s_hat)
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
    isempty(roots) && return (converged=false, choke=true, nu_x=NaN)
    nu_x = prefer_root == :high ? last(roots) : first(roots)
    return (converged=true, choke=false, nu_x=nu_x)
end

"""
    _advance_row!(model, row, aero, row_radius, station_in, station_out, mu, tau, pi, nu_theta, nu_x, m_tip, prefer_root)

Advance one blade row from `station_in` to `station_out` in non-dimensional form.

Key state variables and physical interpretation:
- `nu_x`: axial velocity non-dimensionalized by reference inlet acoustic speed
  (`nu_x = V_x / a0_ref`), stored per station.
- `nu_theta`: tangential (swirl) velocity in the same scaling
  (`nu_theta = V_theta / a0_ref`), stored per station.
- `nu_u`: blade speed at this row/radius in the same scaling
  (`nu_u = U / a0_ref`). For stators, `speed_ratio_to_ref = 0`, so `nu_u = 0`.
- `tau`: non-dimensional total enthalpy / temperature-like state
  (`tau = h_t / h_t_ref` for this ideal-gas closure), stored per station.
- `pi`: non-dimensional total pressure ratio state (`pi = p_t / p_t_ref`),
  stored per station.
- `mu`: mass-flow invariant used by continuity closure. At each stage update we
  solve for `nu_x_out` such that this invariant is preserved.
- `m_tip`: reference rotor-tip speed parameter (`m_tip = omega_ref * r_tip_ref / a0_ref`),
  used to build row-local blade speed via `speed_ratio_to_ref` and `row_radius`.

Behavior:
- Calls `blade_aero` to get turning and loss for the row.
- Solves a scalar continuity residual for `nu_x_out` at the row exit.
- Updates `nu_theta`, `tau`, and `pi`; stator behavior is the `nu_u = 0` special case.
"""
function _advance_row!(
    model::AxialMachineModel,
    row::AxialRow,
    aero::BladeAeroModel,
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
    # nu_u is the blade-relative tangential velocity non-dimensionalized by the reference tip speed.
    # For stators, speed_ratio_to_ref = 0, so nu_u = 0 and the blade-relative flow angle is the same as the absolute flow angle.
    nu_u = row.speed_ratio_to_ref * m_tip * row_radius / model.r_tip_ref
    area_in = station_area(model, station_in)
    area_out = station_area(model, station_out)
    nu_x_in = nu_x[station_in]
    nu_theta_in = nu_theta[station_in]
    tau_in = tau[station_in]
    pi_in = pi[station_in]
    
    # Compute blade aerodynamics from the incoming relative flow and the row's
    # metal inlet/exit angles. Incidence is measured against `theta_metal_in`
    # and deviation is measured from `theta_metal_out`.
    aero_out = blade_aero(
        aero,
        row.theta_metal_in,
        row.theta_metal_out,
        nu_x_in,
        nu_theta_in,
        nu_u,
    )
    stall = (aero_out.stall_margin <= 0) || !aero_out.valid
    diagnostics_base = (
        row_index=station_in,
        station_in=station_in,
        station_out=station_out,
        r_hub=row.r_hub,
        r_tip=row.r_tip,
        row_radius=row_radius,
        row_annulus_area=row_annulus_area(row),
        area_in=area_in,
        area_out=area_out,
        theta_metal_in=row.theta_metal_in,
        theta_metal_out=row.theta_metal_out,
        speed_ratio_to_ref=row.speed_ratio_to_ref,
        nu_u=nu_u,
        nu_x_in=nu_x_in,
        nu_x_out=NaN,
        nu_theta_in=nu_theta_in,
        nu_theta_out=NaN,
        delta_nu_theta=NaN,
        tau_in=tau_in,
        tau_out=NaN,
        delta_tau=NaN,
        pi_in=pi_in,
        pi_out=NaN,
        delta_pi=NaN,
        k_theta_exit=Float64(aero_out.k_theta_exit),
        delta_s_hat=Float64(aero_out.delta_s_hat),
        stall_margin=Float64(aero_out.stall_margin),
        incidence=Float64(aero_out.diagnostics.incidence),
        deviation=Float64(aero_out.diagnostics.deviation),
        theta_in=Float64(aero_out.diagnostics.theta_in),
        theta_out=Float64(aero_out.diagnostics.theta_out),
        valid=Bool(aero_out.valid),
        stall=stall,
        choke=false,
    )
    aero_out.valid || return (converged=false, choke=false, stall=stall, diagnostics=diagnostics_base)

    # Compute nu_x at the row exit from blade aerodynamics and mass flow constraint.
    nu_x_solve = _solve_row_nu_x(
        model,
        station_in,
        station_out,
        mu,
        tau,
        pi,
        nu_theta,
        nu_u,
        Float64(aero_out.k_theta_exit),
        Float64(aero_out.delta_s_hat),
        prefer_root,
    )
    nu_x_solve.converged || return (
        converged=false,
        choke=nu_x_solve.choke,
        stall=stall,
        diagnostics=merge(diagnostics_base, (choke=nu_x_solve.choke,)),
    )
    nu_x_out = nu_x_solve.nu_x

    # Calculate incrase in enthalpy
    gamma_ratio = model.gamma / (model.gamma - 1)
    nu_theta_out = nu_u + nu_x_out * aero_out.k_theta_exit
    tau_out = tau_in + (model.gamma - 1) * nu_u * (nu_theta_out - nu_theta_in)
    tau_out > 0 || return (converged=false, choke=false, stall=stall, diagnostics=diagnostics_base)

    # Calculate increase in pressure ratio
    pi_out = (tau_out / tau_in)^gamma_ratio *
             pi_in * exp(-gamma_ratio * aero_out.delta_s_hat)

    nu_theta[station_out] = nu_theta_out
    tau[station_out] = tau_out
    pi[station_out] = pi_out
    nu_x[station_out] = nu_x_out
    diagnostics = merge(
        diagnostics_base,
        (
            nu_x_out=nu_x_out,
            nu_theta_out=nu_theta_out,
            delta_nu_theta=nu_theta_out - nu_theta_in,
            tau_out=tau_out,
            delta_tau=tau_out - tau_in,
            pi_out=pi_out,
            delta_pi=pi_out - pi_in,
        ),
    )
    return (converged=true, choke=false, stall=stall, diagnostics=diagnostics)
end

"""
    streamtube_solve(model, streamtube_radii, m_tip, nu_x_inlet, nu_theta_inlet; prefer_root=:low)

Run the axial row-marching solve in non-dimensional coordinates.

Coordinate/scaling convention:
- All velocity-like terms are normalized by `a0_ref` (reference inlet acoustic speed).
  - `nu_x = V_x / a0_ref` (axial component)
  - `nu_theta = V_theta / a0_ref` (tangential/swirl component)
  - `nu_u = U / a0_ref` (blade speed, row-dependent)
- Thermodynamic station states are represented by:
  - `tau` (non-dimensional total enthalpy/temperature-like state)
  - `pi` (non-dimensional total pressure ratio state)
- `mu` is a non-dimensional mass-flow invariant derived from continuity.

Inputs:
- `model`: row stack plus gas constants and reference geometry.
- `streamtube_radii[k]`: physical radius used for row `k`; must lie in that row's
  `[r_hub, r_tip]`.
- `m_tip`: reference-speed parameter driving row blade speeds.
- `nu_x_inlet`, `nu_theta_inlet`: inlet station velocity components in normalized units.
- `prefer_root`: selects lower or upper branch when scalar continuity has multiple roots.

Outputs include:
- `PR`: outlet-to-inlet total pressure ratio (`pi_out`)
- `eta`: total-to-total efficiency proxy from `(pi_out, tau_out)`
- full station arrays (`tau`, `pi`, `nu_theta`, `nu_x`) and row diagnostics.
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
    row_data = Vector{NamedTuple}(undef, n_rows)

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
        row_data[k] = row_step.diagnostics
        if !row_step.converged
            choke_row[k] = row_step.choke
            valid_row[k] = false
            break
        end
    end

    for k in 1:n_rows
        if !isassigned(row_data, k)
            row = model.rows[k]
            station_in = k
            station_out = k + 1
            area_in = station_area(model, station_in)
            area_out = station_area(model, station_out)
            row_data[k] = (
                row_index=k,
                station_in=station_in,
                station_out=station_out,
                r_hub=row.r_hub,
                r_tip=row.r_tip,
                row_radius=radii[k],
                row_annulus_area=row_annulus_area(row),
                area_in=area_in,
                area_out=area_out,
                theta_metal_in=row.theta_metal_in,
                theta_metal_out=row.theta_metal_out,
                speed_ratio_to_ref=row.speed_ratio_to_ref,
                nu_u=row.speed_ratio_to_ref * m_tip_f * radii[k] / model.r_tip_ref,
                nu_x_in=nu_x[station_in],
                nu_x_out=nu_x[station_out],
                nu_theta_in=nu_theta[station_in],
                nu_theta_out=nu_theta[station_out],
                delta_nu_theta=NaN,
                tau_in=tau[station_in],
                tau_out=tau[station_out],
                delta_tau=NaN,
                pi_in=pi[station_in],
                pi_out=pi[station_out],
                delta_pi=NaN,
                k_theta_exit=NaN,
                delta_s_hat=NaN,
                stall_margin=NaN,
                incidence=NaN,
                deviation=NaN,
                theta_in=NaN,
                theta_out=NaN,
                valid=false,
                stall=false,
                choke=false,
            )
        end
    end

    station_data = _build_station_data(model, radii, tau, pi, nu_theta, nu_x)

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
            station_data=station_data,
            row_data=row_data,
        )
    end

    tau_out = tau[end]
    pi_out = pi[end]
    eta_tt = if (tau_out > 0) && (pi_out > 0)
        tau_is_out = pi_out^((model.gamma - 1) / model.gamma)
        d_actual = tau_out - 1
        d_is = tau_is_out - 1
        eta = if abs(d_actual) <= 1e-12 || abs(d_is) <= 1e-12
            NaN
        elseif signbit(d_actual) != signbit(d_is)
            NaN
        elseif d_actual > 0
            d_is / d_actual
        else
            d_actual / d_is
        end
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
        station_data=station_data,
        row_data=row_data,
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
    # The phi wrapper uses the model's explicit flow-reference speed/radius
    # to map flow coefficient to inlet axial velocity:
    #   phi_in = nu_x_inlet / |nu_u_ref|  =>  nu_x_inlet = phi_in * |nu_u_ref|
    # where nu_u_ref is the non-dimensional blade speed at `r_flow_ref`.
    nu_u_ref = model.speed_ratio_ref * m_tip_f * model.r_flow_ref / model.r_tip_ref
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

function _nearest_feasible_flow_sample(
    model::AxialMachineModel,
    speed::Float64,
    flow_target::Float64,
    flow_lo::Float64,
    flow_hi::Float64,
    streamtube_radii::AbstractVector{<:Real},
    nu_theta_inlet::Float64,
    prefer_root::Symbol,
    is_feasible::Function;
    n_probe::Int=61,
)
    flow_min = min(flow_lo, flow_hi)
    flow_max = max(flow_lo, flow_hi)
    n_probe >= 3 || error("n_probe must be >= 3")
    probe = collect(range(flow_min, flow_max, length=n_probe))
    push!(probe, flow_target)
    sort!(probe)

    best_dist = Inf
    best_flow = NaN
    best_vals = nothing
    for flow in probe
        vals = streamtube_solve_with_phi(
            model,
            speed,
            flow;
            streamtube_radii=streamtube_radii,
            nu_theta_inlet=nu_theta_inlet,
            prefer_root=prefer_root,
        )
        if is_feasible(vals)
            dist = abs(flow - flow_target)
            if dist < best_dist
                best_dist = dist
                best_flow = flow
                best_vals = vals
            end
        end
    end

    return (
        found=(best_vals !== nothing),
        flow=best_flow,
        vals=best_vals,
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
            if !is_feasible(vals)
                if has_limits
                    lo = Float64(flow_min[i])
                    hi = Float64(flow_max[i])
                    repaired = _nearest_feasible_flow_sample(
                        model,
                        speed,
                        flow,
                        lo,
                        hi,
                        streamtube_radii,
                        Float64(nu_theta_inlet),
                        prefer_root,
                        is_feasible;
                        n_probe=61,
                    )
                    if repaired.found
                        flow = repaired.flow
                        vals = repaired.vals
                    end
                end
            end
            is_feasible(vals) ||
                error("streamtube sampling produced invalid value at speed=$(speed), flow=$(flow), limits=[$(has_limits ? flow_min[i] : NaN), $(has_limits ? flow_max[i] : NaN)]")
            pr_table[i, j] = vals.PR
            eta_table[i, j] = vals.eta
        end
    end

    return (pr_table=pr_table, eta_table=eta_table)
end
