"""
    feasible_flow_limits(model, eos, speed_grid, flow_lo, flow_hi; pt_in, ht_in, boundary_resolution=401, streamtube_radii=meanline_radii(model), Vtheta_inlet=0.0, prefer_root=:low, is_feasible)

Scan physical inlet mass flow at each physical shaft speed and return feasible-flow bounds.
"""
function feasible_flow_limits(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    speed_grid::AbstractVector{<:Real},
    flow_lo::Real,
    flow_hi::Real;
    pt_in::Real,
    ht_in::Real,
    boundary_resolution::Int=401,
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
    is_feasible::Function=(result -> getproperty(result, :valid)),
)
    boundary_resolution >= 2 || error("boundary_resolution must be >= 2")
    flow_hi > flow_lo || error("flow_hi must be > flow_lo")
    length(speed_grid) >= 1 || error("speed_grid must be non-empty")

    flow_probe = collect(range(Float64(flow_lo), Float64(flow_hi), length=boundary_resolution))
    valid_speed_idx = Int[]
    flow_min = Float64[]
    flow_max = Float64[]

    for (i, speed_raw) in pairs(speed_grid)
        speed = Float64(speed_raw)
        feasible_flows = Float64[]
        for flow in flow_probe
            result = streamtube_solve_from_mdot(
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
            if is_feasible(result)
                push!(feasible_flows, flow)
            end
        end
        isempty(feasible_flows) && continue
        push!(valid_speed_idx, i)
        push!(flow_min, first(feasible_flows))
        push!(flow_max, last(feasible_flows))
    end

    return (
        valid_speed_idx=valid_speed_idx,
        flow_min=flow_min,
        flow_max=flow_max,
    )
end
