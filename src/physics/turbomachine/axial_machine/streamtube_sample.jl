function _nearest_feasible_flow_sample(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    omega::Float64,
    flow_target::Float64,
    flow_lo::Float64,
    flow_hi::Float64,
    pt_in::Float64,
    ht_in::Float64,
    streamtube_radii::AbstractVector{<:Real},
    Vtheta_inlet::Float64,
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
        vals = streamtube_solve_from_mdot(
            model,
            eos,
            streamtube_radii,
            pt_in,
            ht_in,
            omega,
            flow,
            Vtheta_inlet;
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
    sample_streamtube_solve(model, eos, speed_grid, flow_grid; pt_in, ht_in, ...)

Sample the dimensional streamtube solver over a Cartesian grid of physical
`(omega, mdot)` coordinates.
"""
function sample_streamtube_solve(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    speed_grid::AbstractVector{<:Real},
    flow_grid::AbstractVector{<:Real};
    pt_in::Real,
    ht_in::Real,
    flow_min::Union{Nothing,AbstractVector{<:Real}}=nothing,
    flow_max::Union{Nothing,AbstractVector{<:Real}}=nothing,
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    Vtheta_inlet::Real=0.0,
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
            vals = streamtube_solve_from_mdot(
                model,
                eos,
                streamtube_radii,
                Float64(pt_in),
                Float64(ht_in),
                speed,
                flow,
                Float64(Vtheta_inlet);
                prefer_root=prefer_root,
            )
            if !is_feasible(vals) && has_limits
                lo = Float64(flow_min[i])
                hi = Float64(flow_max[i])
                repaired = _nearest_feasible_flow_sample(
                    model,
                    eos,
                    speed,
                    flow,
                    lo,
                    hi,
                    Float64(pt_in),
                    Float64(ht_in),
                    streamtube_radii,
                    Float64(Vtheta_inlet),
                    prefer_root,
                    is_feasible;
                    n_probe=61,
                )
                if repaired.found
                    flow = repaired.flow
                    vals = repaired.vals
                end
            end
            is_feasible(vals) ||
                error("streamtube sampling produced invalid value at omega=$(speed), mdot=$(flow), limits=[$(has_limits ? flow_min[i] : NaN), $(has_limits ? flow_max[i] : NaN)]")
            pr_table[i, j] = vals.PR
            eta_table[i, j] = vals.eta
        end
    end

    return (pr_table=pr_table, eta_table=eta_table)
end
