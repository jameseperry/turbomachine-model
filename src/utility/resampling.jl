"""
Generic 1D linear resampling utilities.

These helpers are intended for cases where sampled curve data needs to be
evaluated on a new 1D grid, optionally after sorting and duplicate collapse.
"""

function _deduplicate_mean(
    x::Vector{Float64},
    y::Vector{Float64},
)
    isempty(x) && return (Float64[], Float64[])
    x_out = Float64[x[1]]
    y_out = Float64[y[1]]
    counts = Int[1]
    for i in 2:length(x)
        if isapprox(x[i], x_out[end]; atol=1e-10, rtol=1e-10)
            y_out[end] += y[i]
            counts[end] += 1
        else
            push!(x_out, x[i])
            push!(y_out, y[i])
            push!(counts, 1)
        end
    end
    for i in eachindex(y_out)
        y_out[i] /= counts[i]
    end
    return (x_out, y_out)
end

function _prepare_resample_inputs(
    xgrid::AbstractVector{<:Real},
    values::AbstractVector{<:Real};
    sort_inputs::Bool=false,
    deduplicate::Symbol=:error,
)
    length(xgrid) == length(values) || error("values length must match xgrid length")
    length(xgrid) >= 2 || error("xgrid must have at least 2 points")

    x = Float64.(xgrid)
    y = Float64.(values)

    if sort_inputs
        order = sortperm(x)
        x = x[order]
        y = y[order]
    end

    if deduplicate == :mean
        x, y = _deduplicate_mean(x, y)
    elseif deduplicate != :error
        error("unsupported deduplicate=$(deduplicate)")
    end

    _strictly_increasing_linear(x) || error("xgrid must be strictly increasing after preprocessing")
    return x, y
end

function _resample_extrapolated_value(
    x::Vector{Float64},
    y::Vector{Float64},
    xq::Float64;
    extrapolation::Symbol,
)
    if xq < x[1]
        extrapolation == :nan && return NaN
        extrapolation == :error && error("query point below xgrid domain")
        extrapolation == :clamp || error("unsupported extrapolation=$(extrapolation)")
        return y[1]
    elseif xq > x[end]
        extrapolation == :nan && return NaN
        extrapolation == :error && error("query point above xgrid domain")
        extrapolation == :clamp || error("unsupported extrapolation=$(extrapolation)")
        return y[end]
    end

    idx = searchsortedlast(x, xq)
    idx = clamp(idx, 1, length(x) - 1)
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    if isapprox(x1, x0; atol=1e-12, rtol=1e-12)
        return y0
    end
    t = (xq - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)
end

"""
    resample_linear(xgrid, values, xq; extrapolation=:clamp, sort_inputs=false, deduplicate=:error)

Linearly resample 1D data defined by `(xgrid, values)` at new query locations `xq`.

Keywords:
- `extrapolation`:
  - `:clamp` clamps to edge values
  - `:nan` returns `NaN` outside the input domain
  - `:error` throws for out-of-domain queries
- `sort_inputs=true` sorts the input pairs by `xgrid` before resampling
- `deduplicate=:mean` averages duplicate x-values after sorting
"""
function resample_linear(
    xgrid::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    xq::AbstractVector{<:Real};
    extrapolation::Symbol=:clamp,
    sort_inputs::Bool=false,
    deduplicate::Symbol=:error,
)
    x, y = _prepare_resample_inputs(
        xgrid,
        values;
        sort_inputs=sort_inputs,
        deduplicate=deduplicate,
    )
    return [
        _resample_extrapolated_value(x, y, Float64(q); extrapolation=extrapolation)
        for q in xq
    ]
end

function resample_linear(
    xgrid::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    xq::Real;
    extrapolation::Symbol=:clamp,
    sort_inputs::Bool=false,
    deduplicate::Symbol=:error,
)
    x, y = _prepare_resample_inputs(
        xgrid,
        values;
        sort_inputs=sort_inputs,
        deduplicate=deduplicate,
    )
    return _resample_extrapolated_value(x, y, Float64(xq); extrapolation=extrapolation)
end
