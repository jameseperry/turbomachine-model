function _nearest_feasible_flow_sample(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    omega::Float64,
    Vx_target::Float64,
    Vx_lo::Float64,
    Vx_hi::Float64,
    pt_in::Float64,
    ht_in::Float64,
    streamtube_radii::AbstractVector{<:Real},
    Vtheta_inlet::Float64,
    prefer_root::Symbol,
    is_feasible::Function;
    n_probe::Int=61,
)
    Vx_min = min(Vx_lo, Vx_hi)
    Vx_max = max(Vx_lo, Vx_hi)
    n_probe >= 3 || error("n_probe must be >= 3")
    probe = collect(range(Vx_min, Vx_max, length=n_probe))
    push!(probe, Vx_target)
    sort!(probe)

    best_dist = Inf
    best_Vx = NaN
    best_vals = nothing
    for Vx in probe
        vals = streamtube_solve(
            model,
            eos,
            streamtube_radii,
            pt_in,
            ht_in,
            omega,
            Vx,
            Vtheta_inlet;
            prefer_root=prefer_root,
        )
        if is_feasible(vals)
            dist = abs(Vx - Vx_target)
            if dist < best_dist
                best_dist = dist
                best_Vx = Vx
                best_vals = vals
            end
        end
    end

    return (
        found=(best_vals !== nothing),
        Vx=best_Vx,
        vals=best_vals,
    )
end

"""
    sample_streamtube_solve(model, eos, speed_grid, Vx_grid; pt_in, ht_in, ...)

Sample the dimensional streamtube solver over a Cartesian grid of physical
`(omega, Vx_inlet)` coordinates and return projected map quantities.
"""
function sample_streamtube_solve(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    speed_grid::AbstractVector{<:Real},
    Vx_grid::AbstractVector{<:Real};
    pt_in::Real,
    ht_in::Real,
    Vx_min::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Vx_max::Union{Nothing,AbstractVector{<:Real}}=nothing,
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
    repair_infeasible_samples::Bool=true,
    infeasible_samples_sink::Union{Nothing,AbstractVector}=nothing,
    is_feasible::Function=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta)),
)
    length(speed_grid) >= 1 || error("speed_grid must be non-empty")
    length(Vx_grid) >= 1 || error("Vx_grid must be non-empty")

    has_limits = !isnothing(Vx_min) || !isnothing(Vx_max)
    if has_limits
        isnothing(Vx_min) && error("Vx_min must be provided when using Vx limits")
        isnothing(Vx_max) && error("Vx_max must be provided when using Vx limits")
        length(Vx_min) == length(speed_grid) || error("Vx_min length must match speed_grid")
        length(Vx_max) == length(speed_grid) || error("Vx_max length must match speed_grid")
    end

    pr_table = Matrix{Float64}(undef, length(speed_grid), length(Vx_grid))
    eta_table = Matrix{Float64}(undef, length(speed_grid), length(Vx_grid))
    mdot_table = Matrix{Float64}(undef, length(speed_grid), length(Vx_grid))

    for (i, speed_raw) in pairs(speed_grid)
        speed = Float64(speed_raw)
        for (j, Vx_raw) in pairs(Vx_grid)
            Vx_target = Float64(Vx_raw)
            vals = streamtube_solve(
                model,
                eos,
                streamtube_radii,
                Float64(pt_in),
                Float64(ht_in),
                speed,
                Vx_target,
                Float64(Vtheta_inlet);
                prefer_root=prefer_root,
            )
            if !is_feasible(vals) && !isnothing(infeasible_samples_sink)
                push!(infeasible_samples_sink, (
                    i_speed=i,
                    i_flow=j,
                    omega=speed,
                    Vx_inlet=Vx_target,
                    mdot=vals.stations[1].mdot_station,
                    pt_in=Float64(pt_in),
                    ht_in=Float64(ht_in),
                    Vtheta_inlet=Float64(Vtheta_inlet),
                    Vx_min=has_limits ? Float64(Vx_min[i]) : NaN,
                    Vx_max=has_limits ? Float64(Vx_max[i]) : NaN,
                    pressure_ratio=vals.PR,
                    efficiency=vals.eta,
                    valid=vals.valid,
                    stall=vals.stall,
                    choke=vals.choke,
                    status=vals.status,
                ))
            end
            Vx_used = Vx_target
            if repair_infeasible_samples && !is_feasible(vals) && has_limits
                lo = Float64(Vx_min[i])
                hi = Float64(Vx_max[i])
                repaired = _nearest_feasible_flow_sample(
                    model,
                    eos,
                    speed,
                    Vx_target,
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
                    Vx_used = repaired.Vx
                    vals = repaired.vals
                end
            end
            is_feasible(vals) ||
                error("streamtube sampling produced invalid value at omega=$(speed), Vx_inlet=$(Vx_target), limits=[$(has_limits ? Vx_min[i] : NaN), $(has_limits ? Vx_max[i] : NaN)]")
            pr_table[i, j] = vals.PR
            eta_table[i, j] = vals.eta
            mdot_table[i, j] = vals.stations[1].mdot_station
        end
    end

    return (pr_table=pr_table, eta_table=eta_table, mdot_table=mdot_table)
end
