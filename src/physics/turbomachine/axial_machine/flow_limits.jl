"""
    feasible_flow_limits(model, speed_grid, flow_lo, flow_hi; boundary_resolution=401, prefer_root=:low, is_feasible)

Scan a flow coordinate at each speed and return feasible-flow bounds.

Inputs:
- `model::AxialMachineModel`: axial machine model to evaluate.
- `speed_grid`: speeds to evaluate.
- `flow_lo`, `flow_hi`: scan interval for the flow coordinate.
- `boundary_resolution`: number of probe points on `[flow_lo, flow_hi]`.
- `prefer_root`: root branch selection for `streamtube_solve`.
- `is_feasible(result)`: predicate that marks a solver result as feasible.

Returns:
- `valid_speed_idx`: indices in `speed_grid` with at least one feasible flow.
- `flow_min`: first feasible flow per valid speed (low-flow boundary).
- `flow_max`: last feasible flow per valid speed (high-flow boundary).
"""
function feasible_flow_limits(
    model::AxialMachineModel,
    speed_grid::AbstractVector{<:Real},
    flow_lo::Real,
    flow_hi::Real;
    boundary_resolution::Int=401,
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

    for (i, speed) in pairs(speed_grid)
        feasible_flows = Float64[]
        for flow in flow_probe
            result = streamtube_solve(model, Float64(speed), flow; prefer_root=prefer_root)
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
