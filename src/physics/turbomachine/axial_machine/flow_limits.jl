"""
    feasible_flow_limits(model, speed_grid, flow_lo, flow_hi; boundary_resolution=401, streamtube_radii=meanline_radii(model), nu_theta_inlet=0.0, prefer_root=:low, is_feasible)

Scan inlet flow coefficient at each speed and return feasible-flow bounds.
"""
function feasible_flow_limits(
    model::AxialMachineModel,
    speed_grid::AbstractVector{<:Real},
    flow_lo::Real,
    flow_hi::Real;
    boundary_resolution::Int=401,
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    nu_theta_inlet::Real=0.0,
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
        idx_ref = model.first_rotor_index
        row_ref = model.rows[idx_ref]
        nu_u_ref = row_ref.speed_ratio_to_ref * speed * Float64(streamtube_radii[idx_ref]) / model.r_tip_ref
        abs(nu_u_ref) > 0 || continue
        feasible_flows = Float64[]
        for flow in flow_probe
            nu_x_inlet = flow * abs(nu_u_ref)
            result = streamtube_solve(
                model,
                streamtube_radii,
                speed,
                nu_x_inlet,
                Float64(nu_theta_inlet);
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
