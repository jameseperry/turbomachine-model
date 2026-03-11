"""
    feasible_flow_limits(model, eos, speed_grid, Vx_lo, Vx_hi; pt_in, ht_in, boundary_resolution=401, streamtube_radii=meanline_radii(model), Vtheta_inlet=0.0, prefer_root=:low, is_feasible)

Scan physical inlet axial velocity at each physical shaft speed and return
feasible inlet-velocity bounds, along with the projected mass-flow bounds.
"""
function feasible_flow_limits(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    speed_grid::AbstractVector{<:Real},
    Vx_lo::Real,
    Vx_hi::Real;
    pt_in::Real,
    ht_in::Real,
    boundary_resolution::Int=401,
    streamtube_radii::AbstractVector{<:Real}=meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
    is_feasible::Function=(result -> getproperty(result, :valid)),
)
    boundary_resolution >= 2 || error("boundary_resolution must be >= 2")
    Vx_hi > Vx_lo || error("Vx_hi must be > Vx_lo")
    length(speed_grid) >= 1 || error("speed_grid must be non-empty")

    Vx_probe = collect(range(Float64(Vx_lo), Float64(Vx_hi), length=boundary_resolution))
    valid_speed_idx = Int[]
    Vx_min = Float64[]
    Vx_max = Float64[]
    mdot_min = Float64[]
    mdot_max = Float64[]

    for (i, speed_raw) in pairs(speed_grid)
        speed = Float64(speed_raw)
        feasible_Vx = Float64[]
        feasible_mdot = Float64[]
        for Vx in Vx_probe
            result = streamtube_solve(
                model,
                eos,
                streamtube_radii,
                Float64(pt_in),
                Float64(ht_in),
                speed,
                Vx,
                Float64(Vtheta_inlet);
                prefer_root=prefer_root,
            )
            if is_feasible(result)
                push!(feasible_Vx, Vx)
                push!(feasible_mdot, result.stations[1].mdot_station)
            end
        end
        isempty(feasible_Vx) && continue
        push!(valid_speed_idx, i)
        push!(Vx_min, first(feasible_Vx))
        push!(Vx_max, last(feasible_Vx))
        push!(mdot_min, first(feasible_mdot))
        push!(mdot_max, last(feasible_mdot))
    end

    return (
        valid_speed_idx=valid_speed_idx,
        Vx_min=Vx_min,
        Vx_max=Vx_max,
        mdot_min=mdot_min,
        mdot_max=mdot_max,
    )
end
