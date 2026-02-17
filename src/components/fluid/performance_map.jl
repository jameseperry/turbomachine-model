"""
Map wrapper for compressor/turbine performance data with correction helpers.

Table convention:
- `n_corr_grid` is corrected shaft speed axis.
- `w_corr_grid` is corrected mass-flow axis.
- `pr_table[i, j]` and `eta_table[i, j]` correspond to
  `n_corr_grid[i]`, `w_corr_grid[j]`.
"""
struct PerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    n_corr_grid::Vector{Float64}
    w_corr_grid::Vector{Float64}
    pr_table::Matrix{Float64}
    eta_table::Matrix{Float64}
end

function PerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    n_corr_grid::Vector{<:Real},
    w_corr_grid::Vector{<:Real},
    pr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real},
)
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(n_corr_grid) >= 2 || error("n_corr_grid must have at least 2 points")
    length(w_corr_grid) >= 2 || error("w_corr_grid must have at least 2 points")
    issorted(n_corr_grid) || error("n_corr_grid must be sorted ascending")
    issorted(w_corr_grid) || error("w_corr_grid must be sorted ascending")
    size(pr_table) == (length(n_corr_grid), length(w_corr_grid)) ||
        error("pr_table size must match (length(n_corr_grid), length(w_corr_grid))")
    size(eta_table) == (length(n_corr_grid), length(w_corr_grid)) ||
        error("eta_table size must match (length(n_corr_grid), length(w_corr_grid))")

    return PerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        Float64.(n_corr_grid),
        Float64.(w_corr_grid),
        Float64.(pr_table),
        Float64.(eta_table),
    )
end

corrected_speed(N::Real, Tt_in::Real, Tt_ref::Real) =
    N / sqrt(Tt_in / Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(N::Real, Tt_in::Real, map::PerformanceMap) =
    corrected_speed(N, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::PerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

function _cell_index_and_fraction(grid::Vector{Float64}, x::Float64)
    x_clamped = clamp(x, first(grid), last(grid))
    idx_hi = searchsortedfirst(grid, x_clamped)
    if idx_hi <= 1
        return 1, 0.0
    elseif idx_hi > length(grid)
        return length(grid) - 1, 1.0
    end
    idx_lo = idx_hi - 1
    x0 = grid[idx_lo]
    x1 = grid[idx_hi]
    t = (x_clamped - x0) / (x1 - x0)
    return idx_lo, t
end

function _bilinear(
    x::Float64,
    y::Float64,
    xgrid::Vector{Float64},
    ygrid::Vector{Float64},
    table::Matrix{Float64},
)
    i, tx = _cell_index_and_fraction(xgrid, x)
    j, ty = _cell_index_and_fraction(ygrid, y)

    f00 = table[i, j]
    f10 = table[i + 1, j]
    f01 = table[i, j + 1]
    f11 = table[i + 1, j + 1]

    f0 = f00 * (1 - tx) + f10 * tx
    f1 = f01 * (1 - tx) + f11 * tx
    return f0 * (1 - ty) + f1 * ty
end

function map_pr_eta(map::PerformanceMap, N_corr::Real, W_corr::Real)
    n = Float64(N_corr)
    w = Float64(W_corr)
    PR = _bilinear(n, w, map.n_corr_grid, map.w_corr_grid, map.pr_table)
    eta = _bilinear(n, w, map.n_corr_grid, map.w_corr_grid, map.eta_table)
    return (PR=PR, eta=eta)
end

function map_pr_eta_from_inlet(
    map::PerformanceMap,
    N::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    N_corr = corrected_speed(N, Tt_in, map)
    W_corr = corrected_flow(mdot, Tt_in, Pt_in, map)
    vals = map_pr_eta(map, N_corr, W_corr)
    return (N_corr=N_corr, W_corr=W_corr, PR=vals.PR, eta=vals.eta)
end
