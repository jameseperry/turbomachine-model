"""
Operating-point sweep built on top of the branch-complete point solver.
"""

using ..Fluids
using ...Utility: track_branches

function _resolve_pt_out_grid(
    pt_in::Float64;
    pt_out_grid::Union{Nothing,AbstractVector}=nothing,
    pr_grid::Union{Nothing,AbstractVector}=nothing,
)
    if isnothing(pt_out_grid) == isnothing(pr_grid)
        error("exactly one of pt_out_grid or pr_grid must be provided")
    end
    if !isnothing(pt_out_grid)
        return Float64.(pt_out_grid)
    end
    return pt_in .* Float64.(pr_grid)
end

function _candidate_branch_coordinates(point_result)
    return [Float64(candidate.branch_coordinate) for candidate in point_result.candidates]
end

function _neighbor_continuation_hints(point_results, i::Int, j::Int)
    hints = Float64[]
    if i > 1
        above = point_results[i - 1, j]
        isnothing(above) || append!(hints, _candidate_branch_coordinates(above))
    end
    if j > 1
        left = point_results[i, j - 1]
        isnothing(left) || append!(hints, _candidate_branch_coordinates(left))
    end
    sort!(hints)
    return unique(hints)
end

function _neighbor_reference_coordinate(
    selected_branch_coordinate::Matrix{Float64},
    i::Int,
    j::Int,
    initial_branch_coordinate::Union{Nothing,Real},
)
    if j > 1 && isfinite(selected_branch_coordinate[i, j - 1])
        return selected_branch_coordinate[i, j - 1]
    end
    if i > 1 && isfinite(selected_branch_coordinate[i - 1, j])
        return selected_branch_coordinate[i - 1, j]
    end
    return isnothing(initial_branch_coordinate) ? nothing : Float64(initial_branch_coordinate)
end

function _select_branch_for_sweep(
    point_result;
    branch_selection::Symbol,
    reference_coordinate::Union{Nothing,Real}=nothing,
)
    isempty(point_result.candidates) && return nothing
    if branch_selection == :nearest && isnothing(reference_coordinate)
        return select_operating_point_branch(point_result; policy=:low)
    end
    if branch_selection in (:low, :high)
        return select_operating_point_branch(point_result; policy=branch_selection)
    elseif branch_selection == :nearest
        return select_operating_point_branch(
            point_result;
            policy=:nearest,
            reference_coordinate=reference_coordinate,
        )
    else
        error("unsupported branch_selection=$(branch_selection)")
    end
end

function _track_sweep_branches(
    omega_grid::Vector{Float64},
    point_results,
    pt_out_grid::Vector{Float64};
    max_match_cost::Float64=Inf,
)
    n_omega = length(omega_grid)
    n_pt_out = length(pt_out_grid)
    tracking_by_pt_out = Vector{NamedTuple}(undef, n_pt_out)
    branch_ids = Matrix{Vector{Int}}(undef, n_omega, n_pt_out)
    for i in 1:n_omega, j in 1:n_pt_out
        branch_ids[i, j] = Int[]
    end

    next_global_branch_id = 1
    for j in 1:n_pt_out
        roots_by_condition = [point_results[i, j].candidates for i in 1:n_omega]
        tracking = track_branches(
            omega_grid,
            roots_by_condition;
            distance=(a, b) -> abs(a.branch_coordinate - b.branch_coordinate),
            max_match_cost=max_match_cost,
        )

        remapped_branches = Dict{Int, Vector{NamedTuple}}()
        id_map = Dict{Int, Int}()
        for local_id in tracking.branch_ids
            global_id = next_global_branch_id
            next_global_branch_id += 1
            id_map[local_id] = global_id
            remapped_branches[global_id] = tracking.branches[local_id]
        end

        remapped_assignments = Vector{Vector{Int}}(undef, length(tracking.assignments))
        for i in eachindex(tracking.assignments)
            remapped_assignments[i] = [id_map[local_id] for local_id in tracking.assignments[i]]
            branch_ids[i, j] = remapped_assignments[i]
        end

        remapped_unmatched = [
            merge(
                point,
                (branch_id=(point.branch_id == 0 ? 0 : id_map[point.branch_id]),),
            ) for point in tracking.unmatched
        ]

        tracking_by_pt_out[j] = (
            pt_out=pt_out_grid[j],
            assignments=remapped_assignments,
            branches=remapped_branches,
            branch_ids=sort!(collect(keys(remapped_branches))),
            unmatched=remapped_unmatched,
        )
    end

    return (tracking_by_pt_out=tracking_by_pt_out, branch_ids=branch_ids)
end

"""
    solve_operating_sweep(map, eos; ...)

Solve a grid of operating points by repeatedly calling `solve_operating_point(...)`.
The point solver returns all admissible branches; the sweep layer either tracks all
branches or selects one branch per sweep point.
"""
function solve_operating_sweep(
    map,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    omega_grid::AbstractVector,
    pt_out_grid::Union{Nothing,AbstractVector}=nothing,
    pr_grid::Union{Nothing,AbstractVector}=nothing,
    branch_selection::Symbol=:nearest,
    initial_branch_coordinate::Union{Nothing,Real}=nothing,
    track_all_branches::Bool=false,
    want_diagnostics::Bool=true,
    options::NamedTuple=NamedTuple(),
    max_match_cost::Float64=Inf,
)
    branch_selection in (:low, :high, :nearest) || error("branch_selection must be :low, :high, or :nearest")
    pt_in_f = Float64(pt_in)
    ht_in_f = Float64(ht_in)
    omega_grid_f = Float64.(omega_grid)
    isempty(omega_grid_f) && error("omega_grid must not be empty")
    pt_out_grid_f = _resolve_pt_out_grid(pt_in_f; pt_out_grid=pt_out_grid, pr_grid=pr_grid)
    isempty(pt_out_grid_f) && error("pt_out_grid/pr_grid must not be empty")

    n_omega = length(omega_grid_f)
    n_pt_out = length(pt_out_grid_f)
    point_results = Matrix{Any}(undef, n_omega, n_pt_out)

    for i in 1:n_omega
        for j in 1:n_pt_out
            continuation_hints = _neighbor_continuation_hints(point_results, i, j)
            point_results[i, j] = solve_operating_point(
                map,
                eos;
                pt_in=pt_in_f,
                ht_in=ht_in_f,
                pt_out=pt_out_grid_f[j],
                omega=omega_grid_f[i],
                continuation_hints=continuation_hints,
                want_diagnostics=want_diagnostics,
                options=options,
            )
        end
    end

    diagnostics = Matrix{NamedTuple}(undef, n_omega, n_pt_out)
    for i in 1:n_omega, j in 1:n_pt_out
        diagnostics[i, j] = point_results[i, j].diagnostics
    end

    if track_all_branches
        tracking = _track_sweep_branches(
            omega_grid_f,
            point_results,
            pt_out_grid_f;
            max_match_cost=max_match_cost,
        )
        rows = NamedTuple[]
        for i in 1:n_omega
            for j in 1:n_pt_out
                result = point_results[i, j]
                if isempty(result.candidates)
                    push!(rows, (
                        i_omega=i,
                        i_pt_out=j,
                        omega=omega_grid_f[i],
                        pt_out=pt_out_grid_f[j],
                        branch_id=0,
                        branch_coordinate=NaN,
                        mdot=NaN,
                        ht_out=NaN,
                        tau=NaN,
                        PR=NaN,
                        eta=NaN,
                        power=NaN,
                        converged=false,
                        machine_payload=NamedTuple(),
                    ))
                    continue
                end
                branch_ids = tracking.branch_ids[i, j]
                for (k, candidate) in enumerate(result.candidates)
                    push!(rows, (
                        i_omega=i,
                        i_pt_out=j,
                        omega=omega_grid_f[i],
                        pt_out=pt_out_grid_f[j],
                        branch_id=branch_ids[k],
                        branch_coordinate=candidate.branch_coordinate,
                        mdot=candidate.mdot,
                        ht_out=candidate.ht_out,
                        tau=candidate.tau,
                        PR=candidate.pressure_ratio,
                        eta=candidate.efficiency,
                        power=candidate.power,
                        converged=true,
                        machine_payload=candidate.machine_payload,
                    ))
                end
            end
        end
        return (
            mode=:all,
            omega_grid=omega_grid_f,
            pt_out_grid=pt_out_grid_f,
            point_results=point_results,
            rows=rows,
            tracking=tracking.tracking_by_pt_out,
            diagnostics=diagnostics,
        )
    end

    mdot = fill(NaN, n_omega, n_pt_out)
    ht_out = fill(NaN, n_omega, n_pt_out)
    tau = fill(NaN, n_omega, n_pt_out)
    PR = fill(NaN, n_omega, n_pt_out)
    eta = fill(NaN, n_omega, n_pt_out)
    power = fill(NaN, n_omega, n_pt_out)
    branch_coordinate = fill(NaN, n_omega, n_pt_out)
    branch_index = fill(0, n_omega, n_pt_out)
    converged = falses(n_omega, n_pt_out)

    for i in 1:n_omega
        for j in 1:n_pt_out
            result = point_results[i, j]
            ref = _neighbor_reference_coordinate(branch_coordinate, i, j, initial_branch_coordinate)
            selected = _select_branch_for_sweep(
                result;
                branch_selection=branch_selection,
                reference_coordinate=ref,
            )
            isnothing(selected) && continue
            mdot[i, j] = selected.mdot
            ht_out[i, j] = selected.ht_out
            tau[i, j] = selected.tau
            PR[i, j] = selected.PR
            eta[i, j] = selected.eta
            power[i, j] = selected.power
            branch_coordinate[i, j] = selected.branch_coordinate
            branch_index[i, j] = selected.branch_index
            converged[i, j] = true
        end
    end

    return (
        mode=:single,
        omega_grid=omega_grid_f,
        pt_out_grid=pt_out_grid_f,
        mdot=mdot,
        ht_out=ht_out,
        tau=tau,
        PR=PR,
        eta=eta,
        power=power,
        branch_coordinate=branch_coordinate,
        branch_index=branch_index,
        converged=converged,
        diagnostics=diagnostics,
        point_results=point_results,
    )
end
