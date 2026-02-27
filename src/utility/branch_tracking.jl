"""
Generic branch tracking for sets of roots across ordered operating conditions.
"""

function _greedy_min_cost_pairs(costs::Matrix{Float64}, max_cost::Float64)
    n_rows, n_cols = size(costs)
    used_rows = falses(n_rows)
    used_cols = falses(n_cols)
    pairs = NamedTuple[]

    while true
        best_cost = Inf
        best_row = 0
        best_col = 0

        for i in 1:n_rows
            used_rows[i] && continue
            for j in 1:n_cols
                used_cols[j] && continue
                c = costs[i, j]
                if c < best_cost
                    best_cost = c
                    best_row = i
                    best_col = j
                end
            end
        end

        if !(isfinite(best_cost) && best_cost <= max_cost) || best_row == 0
            break
        end

        used_rows[best_row] = true
        used_cols[best_col] = true
        push!(pairs, (row=best_row, col=best_col, cost=best_cost))
    end

    return pairs, used_rows, used_cols
end

"""
    track_branches(x, roots_by_condition; distance, max_match_cost, allow_birth=true, allow_death=true)

Track persistent branches across an ordered list of operating conditions.

Arguments:
- `x`: ordered condition values (for example, shaft speed sweep points).
- `roots_by_condition[i]`: root values found at `x[i]`.
- `distance(a, b)`: user-supplied nonnegative cost for matching two roots.
- `max_match_cost`: maximum allowed match cost.

Returns a named tuple:
- `assignments`: `assignments[i][j]` is the persistent branch id for the j-th root at i-th condition.
- `branches`: dictionary `branch_id => points`, where each point is
  `(x, condition_idx, local_idx, root, match_cost)`.
- `branch_ids`: sorted list of branch ids.
- `unmatched`: diagnostics for dropped/unassigned roots.
"""
function track_branches(
    x::AbstractVector{<:Real},
    roots_by_condition::AbstractVector{<:AbstractVector};
    distance::Function,
    max_match_cost::Real,
    allow_birth::Bool=true,
    allow_death::Bool=true,
)
    length(x) == length(roots_by_condition) || error("x and roots_by_condition must have the same length")
    max_match_cost_f = Float64(max_match_cost)
    max_match_cost_f >= 0 || error("max_match_cost must be >= 0")

    n = length(x)
    assignments = Vector{Vector{Int}}(undef, n)
    branches = Dict{Int, Vector{NamedTuple}}()
    active_root_by_branch = Dict{Int, Any}()
    unmatched = NamedTuple[]
    next_branch_id = 1

    for i in 1:n
        xi = Float64(x[i])
        roots_i = roots_by_condition[i]
        n_roots = length(roots_i)
        assigned_ids = fill(0, n_roots)

        if i == 1
            for j in 1:n_roots
                root = roots_i[j]
                bid = next_branch_id
                next_branch_id += 1
                assigned_ids[j] = bid
                branches[bid] = [(x=xi, condition_idx=i, local_idx=j, root=root, match_cost=NaN)]
                active_root_by_branch[bid] = root
            end
            assignments[i] = assigned_ids
            continue
        end

        active_ids = sort!(collect(keys(active_root_by_branch)))
        if !isempty(active_ids) && n_roots > 0
            costs = fill(Inf, length(active_ids), n_roots)
            for (row, bid) in enumerate(active_ids)
                prev_root = active_root_by_branch[bid]
                for col in 1:n_roots
                    c = Float64(distance(prev_root, roots_i[col]))
                    costs[row, col] = c
                end
            end

            pairs, used_rows, used_cols = _greedy_min_cost_pairs(costs, max_match_cost_f)

            for pair in pairs
                bid = active_ids[pair.row]
                local_idx = pair.col
                root = roots_i[local_idx]
                assigned_ids[local_idx] = bid
                push!(
                    branches[bid],
                    (
                        x=xi,
                        condition_idx=i,
                        local_idx=local_idx,
                        root=root,
                        match_cost=pair.cost,
                    ),
                )
                active_root_by_branch[bid] = root
            end

            for (row, bid) in enumerate(active_ids)
                if !used_rows[row] && allow_death
                    push!(
                        unmatched,
                        (
                            kind=:death,
                            condition_idx=i,
                            branch_id=bid,
                            local_idx=0,
                            x=xi,
                            root=active_root_by_branch[bid],
                        ),
                    )
                    delete!(active_root_by_branch, bid)
                end
            end

            for col in 1:n_roots
                if used_cols[col]
                    continue
                end
                root = roots_i[col]
                if allow_birth
                    bid = next_branch_id
                    next_branch_id += 1
                    assigned_ids[col] = bid
                    branches[bid] = [
                        (
                            x=xi,
                            condition_idx=i,
                            local_idx=col,
                            root=root,
                            match_cost=NaN,
                        ),
                    ]
                    active_root_by_branch[bid] = root
                else
                    push!(
                        unmatched,
                        (
                            kind=:unassigned,
                            condition_idx=i,
                            branch_id=0,
                            local_idx=col,
                            x=xi,
                            root=root,
                        ),
                    )
                end
            end
        elseif n_roots > 0
            for col in 1:n_roots
                root = roots_i[col]
                if allow_birth
                    bid = next_branch_id
                    next_branch_id += 1
                    assigned_ids[col] = bid
                    branches[bid] = [
                        (
                            x=xi,
                            condition_idx=i,
                            local_idx=col,
                            root=root,
                            match_cost=NaN,
                        ),
                    ]
                    active_root_by_branch[bid] = root
                else
                    push!(
                        unmatched,
                        (
                            kind=:unassigned,
                            condition_idx=i,
                            branch_id=0,
                            local_idx=col,
                            x=xi,
                            root=root,
                        ),
                    )
                end
            end
        else
            if allow_death
                for bid in active_ids
                    push!(
                        unmatched,
                        (
                            kind=:death,
                            condition_idx=i,
                            branch_id=bid,
                            local_idx=0,
                            x=xi,
                            root=active_root_by_branch[bid],
                        ),
                    )
                    delete!(active_root_by_branch, bid)
                end
            end
        end

        assignments[i] = assigned_ids
    end

    return (
        assignments=assignments,
        branches=branches,
        branch_ids=sort!(collect(keys(branches))),
        unmatched=unmatched,
    )
end
