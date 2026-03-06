"""
Generic 1D linear interpolation map on a strictly increasing grid.

Fields:
- `xgrid`: x-axis coordinates (ascending).
- `values`: function values at each x-grid point.
"""
struct LinearMap
    xgrid::Vector{Float64}
    values::Vector{Float64}
    function LinearMap(
        xgrid::Vector{Float64},
        values::Vector{Float64},
    )
        length(xgrid) >= 2 || error("xgrid must have at least 2 points")
        length(values) == length(xgrid) || error("values length must match xgrid length")
        _strictly_increasing_linear(xgrid) || error("xgrid must be strictly increasing")
        return new(xgrid, values)
    end
end

@inline function _primal_value_linear(x::Real)
    return hasfield(typeof(x), :value) ? getfield(x, :value) : x
end

_strictly_increasing_linear(v::AbstractVector{<:Real}) =
    all(v[i + 1] > v[i] for i in 1:(length(v) - 1))

function LinearMap(
    xgrid::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
)
    return LinearMap(Float64.(xgrid), Float64.(values))
end

@inline function _linear_index_fraction(
    xgrid::AbstractVector{<:Real},
    x::Real,
)
    x_primal = _primal_value_linear(x)
    if x_primal <= first(xgrid)
        return 1, zero(x)
    elseif x_primal >= last(xgrid)
        return length(xgrid) - 1, one(x)
    end
    idx_hi = searchsortedfirst(xgrid, x_primal)
    idx_lo = idx_hi - 1
    x0 = xgrid[idx_lo]
    x1 = xgrid[idx_hi]
    t = (x - x0) / (x1 - x0)
    return idx_lo, t
end

"""
Evaluate `map` by 1D linear interpolation with clamped extrapolation to edge
segments.
"""
function linear_evaluate(map::LinearMap, x::Real)
    i, t = _linear_index_fraction(map.xgrid, x)
    y0 = map.values[i]
    y1 = map.values[i + 1]
    return y0 * (1 - t) + y1 * t
end

"""
Convenience vector API for clamped 1D linear interpolation.

Validates the input grid shape and monotonicity before evaluating.
"""
function linear_evaluate(
    xgrid::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    x::Real,
)
    length(xgrid) >= 2 || error("xgrid must have at least 2 points")
    length(values) == length(xgrid) || error("values length must match xgrid length")
    _strictly_increasing_linear(xgrid) || error("xgrid must be strictly increasing")
    i, t = _linear_index_fraction(xgrid, x)
    y0 = values[i]
    y1 = values[i + 1]
    return y0 * (1 - t) + y1 * t
end
