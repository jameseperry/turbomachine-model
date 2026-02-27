"""
Generic scalar root-finding and feasibility-backoff utilities.
"""

function _dedupe_sorted_values(values::Vector{Float64}; atol::Float64)
    isempty(values) && return values
    sort!(values)
    out = Float64[values[1]]
    for i in 2:length(values)
        if abs(values[i] - out[end]) > atol
            push!(out, values[i])
        end
    end
    return out
end

function _bisect_zero(
    f::Function,
    a::Float64,
    b::Float64;
    tol::Float64=1e-8,
    max_iters::Int=60,
)
    fa = f(a)
    fb = f(b)
    fa * fb <= 0 || error("invalid bisection bracket")

    lo = a
    hi = b
    flo = fa
    mid = 0.5 * (lo + hi)
    for _ in 1:max_iters
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) <= tol || abs(hi - lo) <= 1e-12
            break
        end
        if flo * fmid <= 0
            hi = mid
        else
            lo = mid
            flo = fmid
        end
    end
    return mid
end

"""
    bracket_bisect_roots(f, (x_lo, x_hi); n_scan=401, root_tol=1e-8, prior_roots=[], ...)

Generic scalar multi-root finder based on sign-change bracketing and bisection.
`f(x)` is the residual whose zeros are sought. `prior_roots` are optional
continuation hints that get injected into the scan grid (with small neighborhoods)
to improve root tracking between nearby conditions.
"""
function bracket_bisect_roots(
    f::Function,
    x_range::Tuple{Float64,Float64};
    n_scan::Int=401,
    root_tol::Float64=1e-8,
    prior_roots::AbstractVector{<:Real}=Float64[],
    continuation_band_fraction::Float64=0.02,
    max_bisect_iters::Int=60,
)
    n_scan >= 5 || error("n_scan must be >= 5")
    x_lo, x_hi = x_range
    x_hi > x_lo || error("x_range must satisfy x_hi > x_lo")

    grid = collect(range(x_lo, x_hi, length=n_scan))
    band = continuation_band_fraction * (x_hi - x_lo)
    for r_raw in prior_roots
        r = Float64(r_raw)
        if x_lo <= r <= x_hi
            push!(grid, r)
            push!(grid, clamp(r - band, x_lo, x_hi))
            push!(grid, clamp(r + band, x_lo, x_hi))
        end
    end
    sort!(grid)
    grid = unique(grid)

    f_vals = Vector{Float64}(undef, length(grid))
    for i in eachindex(grid)
        f_vals[i] = f(grid[i])
    end

    roots = Float64[]
    for i in 1:(length(grid) - 1)
        x1 = grid[i]
        x2 = grid[i + 1]
        f1 = f_vals[i]
        f2 = f_vals[i + 1]
        if abs(f1) <= root_tol
            push!(roots, x1)
            continue
        end
        if f1 * f2 < 0.0
            push!(roots, _bisect_zero(f, x1, x2; tol=root_tol, max_iters=max_bisect_iters))
        end
    end
    if abs(f_vals[end]) <= root_tol
        push!(roots, grid[end])
    end

    dedupe_tol = max((x_hi - x_lo) / 1e6, 1e-10)
    return _dedupe_sorted_values(roots; atol=dedupe_tol)
end

"""
    feasibility_backoff(evaluate_at, target_value; ...)

Generic 1D feasibility backoff for monotone degradation searches.

Arguments:
- `evaluate_at(value)` returns any result object for that value.
- `is_feasible(result)` predicate determines feasibility.

Behavior:
1. Try `target_value` (clamped by `max_value`).
2. If infeasible and enabled, do coarse scan down to `min_value`.
3. Refine highest feasible value with bisection.
"""
function feasibility_backoff(
    evaluate_at::Function,
    target_value::Float64;
    enabled::Bool=true,
    min_value::Float64,
    max_value::Float64=target_value,
    value_tol::Float64=1e-6,
    max_iters::Int=24,
    n_probe::Int=33,
    is_feasible::Function,
)
    value_tol > 0 || error("value_tol must be > 0")
    max_iters >= 1 || error("max_iters must be >= 1")
    n_probe >= 3 || error("n_probe must be >= 3")

    lo = min_value
    hi = min(max_value, target_value)
    lo <= hi || return (converged=false, value=NaN, result=nothing, used_backoff=false)

    result_hi = evaluate_at(hi)
    if is_feasible(result_hi)
        return (
            converged=true,
            value=hi,
            result=result_hi,
            used_backoff=(hi < target_value),
        )
    end

    enabled || return (converged=false, value=NaN, result=nothing, used_backoff=false)

    probe_vals = collect(range(hi, lo, length=n_probe))
    found_feasible = false
    feasible_val = NaN
    feasible_result = nothing
    infeasible_above_val = hi

    for v in probe_vals[2:end]
        r = evaluate_at(v)
        if is_feasible(r)
            found_feasible = true
            feasible_val = v
            feasible_result = r
            break
        end
        infeasible_above_val = v
    end

    found_feasible || return (converged=false, value=NaN, result=nothing, used_backoff=false)

    lo_feas = feasible_val
    hi_infeas = infeasible_above_val
    best_val = lo_feas
    best_result = feasible_result

    for _ in 1:max_iters
        if (hi_infeas - lo_feas) <= value_tol
            break
        end
        mid = 0.5 * (lo_feas + hi_infeas)
        r_mid = evaluate_at(mid)
        if is_feasible(r_mid)
            lo_feas = mid
            best_val = mid
            best_result = r_mid
        else
            hi_infeas = mid
        end
    end

    return (converged=true, value=best_val, result=best_result, used_backoff=true)
end
