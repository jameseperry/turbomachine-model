module Utility

export BilinearMap, bilinear_evaluate

"""
Generic bilinear interpolation map on a rectilinear grid.

Fields:
- `xgrid`: x-axis coordinates (ascending).
- `ygrid`: y-axis coordinates (ascending).
- `table`: function values with shape `(length(xgrid), length(ygrid))`.
"""
struct BilinearMap
    xgrid::Vector{Float64}
    ygrid::Vector{Float64}
    table::Matrix{Float64}
end

function BilinearMap(
    xgrid::Vector{<:Real},
    ygrid::Vector{<:Real},
    table::Matrix{<:Real},
)
    length(xgrid) >= 2 || error("xgrid must have at least 2 points")
    length(ygrid) >= 2 || error("ygrid must have at least 2 points")
    issorted(xgrid) || error("xgrid must be sorted ascending")
    issorted(ygrid) || error("ygrid must be sorted ascending")
    size(table) == (length(xgrid), length(ygrid)) ||
        error("table size must match (length(xgrid), length(ygrid))")

    return BilinearMap(Float64.(xgrid), Float64.(ygrid), Float64.(table))
end

@inline _primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x

function _cell_index_and_fraction(grid::AbstractVector{<:Real}, x::Real)
    x_primal = _primal_value(x)
    if x_primal <= first(grid)
        return 1, zero(x)
    elseif x_primal >= last(grid)
        return length(grid) - 1, one(x)
    end
    idx_hi = searchsortedfirst(grid, x_primal)
    idx_lo = idx_hi - 1
    x0 = grid[idx_lo]
    x1 = grid[idx_hi]
    t = (x - x0) / (x1 - x0)
    return idx_lo, t
end

"""
Evaluate `map` by bilinear interpolation with linear extrapolation clamped
to edge cells.
"""
function bilinear_evaluate(map::BilinearMap, x::Real, y::Real)
    i, tx = _cell_index_and_fraction(map.xgrid, x)
    j, ty = _cell_index_and_fraction(map.ygrid, y)

    f00 = map.table[i, j]
    f10 = map.table[i + 1, j]
    f01 = map.table[i, j + 1]
    f11 = map.table[i + 1, j + 1]

    f0 = f00 * (1 - tx) + f10 * tx
    f1 = f01 * (1 - tx) + f11 * tx
    return f0 * (1 - ty) + f1 * ty
end

end # module Utility
