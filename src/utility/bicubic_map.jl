"""
Generic bicubic Hermite interpolation map on a rectilinear grid.

Node derivatives are estimated with monotone-limited 1D slopes to reduce
overshoot relative to unconstrained cubic splines.
"""
struct BicubicMap <: AbstractTableMap
    xgrid::Vector{Float64}
    ygrid::Vector{Float64}
    table::Matrix{Float64}
    fx::Matrix{Float64}
    fy::Matrix{Float64}
    fxy::Matrix{Float64}
end

function BicubicMap(
    xgrid::Vector{<:Real},
    ygrid::Vector{<:Real},
    table::Matrix{<:Real},
)
    _validate_grid_table(xgrid, ygrid, table)

    xgrid_f = Float64.(xgrid)
    ygrid_f = Float64.(ygrid)
    table_f = Float64.(table)
    fx, fy, fxy = _bicubic_partials(xgrid_f, ygrid_f, table_f)
    return BicubicMap(xgrid_f, ygrid_f, table_f, fx, fy, fxy)
end

function _monotone_slopes_1d(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
)::Vector{Float64}
    n = length(x)
    n == length(y) || error("x/y size mismatch")
    n >= 2 || error("x/y must have at least 2 points")

    x_f = Float64.(x)
    y_f = Float64.(y)
    h = [x_f[i + 1] - x_f[i] for i in 1:(n - 1)]
    δ = [(y_f[i + 1] - y_f[i]) / h[i] for i in 1:(n - 1)]

    if n == 2
        return [δ[1], δ[1]]
    end

    m = zeros(Float64, n)

    for i in 2:(n - 1)
        if δ[i - 1] == 0.0 || δ[i] == 0.0 || sign(δ[i - 1]) != sign(δ[i])
            m[i] = 0.0
        else
            w1 = 2.0 * h[i] + h[i - 1]
            w2 = h[i] + 2.0 * h[i - 1]
            m[i] = (w1 + w2) / ((w1 / δ[i - 1]) + (w2 / δ[i]))
        end
    end

    m[1] = ((2.0 * h[1] + h[2]) * δ[1] - h[1] * δ[2]) / (h[1] + h[2])
    if sign(m[1]) != sign(δ[1])
        m[1] = 0.0
    elseif sign(δ[1]) != sign(δ[2]) && abs(m[1]) > 3.0 * abs(δ[1])
        m[1] = 3.0 * δ[1]
    end

    m[n] = ((2.0 * h[n - 1] + h[n - 2]) * δ[n - 1] - h[n - 1] * δ[n - 2]) / (h[n - 1] + h[n - 2])
    if sign(m[n]) != sign(δ[n - 1])
        m[n] = 0.0
    elseif sign(δ[n - 1]) != sign(δ[n - 2]) && abs(m[n]) > 3.0 * abs(δ[n - 1])
        m[n] = 3.0 * δ[n - 1]
    end

    return m
end

function _bicubic_partials(
    xgrid::Vector{Float64},
    ygrid::Vector{Float64},
    table::Matrix{Float64},
)
    nx, ny = size(table)
    fx = zeros(Float64, nx, ny)
    fy = zeros(Float64, nx, ny)
    fxy_a = zeros(Float64, nx, ny)
    fxy_b = zeros(Float64, nx, ny)

    for j in 1:ny
        fx[:, j] = _monotone_slopes_1d(xgrid, view(table, :, j))
    end
    for i in 1:nx
        fy[i, :] = _monotone_slopes_1d(ygrid, vec(view(table, i, :)))
    end
    for i in 1:nx
        fxy_a[i, :] = _monotone_slopes_1d(ygrid, vec(view(fx, i, :)))
    end
    for j in 1:ny
        fxy_b[:, j] = _monotone_slopes_1d(xgrid, view(fy, :, j))
    end

    fxy = 0.5 .* (fxy_a .+ fxy_b)
    return fx, fy, fxy
end

@inline _h00(t) = 2 * t^3 - 3 * t^2 + 1
@inline _h10(t) = t^3 - 2 * t^2 + t
@inline _h01(t) = -2 * t^3 + 3 * t^2
@inline _h11(t) = t^3 - t^2

@inline function _cubic_hermite(t, v0, v1, m0, m1)
    return _h00(t) * v0 + _h10(t) * m0 + _h01(t) * v1 + _h11(t) * m1
end

function bicubic_evaluate(map::BicubicMap, x::Real, y::Real)
    i, tx = _cell_index_and_fraction(map.xgrid, x)
    j, ty = _cell_index_and_fraction(map.ygrid, y)
    dx = map.xgrid[i + 1] - map.xgrid[i]
    dy = map.ygrid[j + 1] - map.ygrid[j]

    f00 = map.table[i, j]
    f10 = map.table[i + 1, j]
    f01 = map.table[i, j + 1]
    f11 = map.table[i + 1, j + 1]

    fx00 = map.fx[i, j]
    fx10 = map.fx[i + 1, j]
    fx01 = map.fx[i, j + 1]
    fx11 = map.fx[i + 1, j + 1]

    fy00 = map.fy[i, j]
    fy10 = map.fy[i + 1, j]
    fy01 = map.fy[i, j + 1]
    fy11 = map.fy[i + 1, j + 1]

    fxy00 = map.fxy[i, j]
    fxy10 = map.fxy[i + 1, j]
    fxy01 = map.fxy[i, j + 1]
    fxy11 = map.fxy[i + 1, j + 1]

    g0 = _cubic_hermite(tx, f00, f10, fx00 * dx, fx10 * dx)
    g1 = _cubic_hermite(tx, f01, f11, fx01 * dx, fx11 * dx)
    gy0 = _cubic_hermite(tx, fy00, fy10, fxy00 * dx, fxy10 * dx)
    gy1 = _cubic_hermite(tx, fy01, fy11, fxy01 * dx, fxy11 * dx)

    return _cubic_hermite(ty, g0, g1, gy0 * dy, gy1 * dy)
end

table_evaluate(map::BicubicMap, x::Real, y::Real) = bicubic_evaluate(map, x, y)
table_xgrid(map::BicubicMap) = map.xgrid
table_ygrid(map::BicubicMap) = map.ygrid
table_values(map::BicubicMap) = map.table
table_interpolation(::BicubicMap) = :bicubic
