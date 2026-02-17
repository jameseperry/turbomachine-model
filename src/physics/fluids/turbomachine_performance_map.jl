"""
Map wrapper for compressor/turbine performance data with correction helpers.

Table convention:
- `n_corr_grid` is corrected shaft speed axis.
- `w_corr_grid` is corrected mass-flow axis.
- `pr_table[i, j]` and `eta_table[i, j]` correspond to
  `n_corr_grid[i]`, `w_corr_grid[j]`.

Parameters:
- `Tt_ref`: reference total temperature used for corrected-value normalization.
- `Pt_ref`: reference total pressure used for corrected-value normalization.
- `n_corr_grid`: corrected speed grid (ascending).
- `w_corr_grid`: corrected mass-flow grid (ascending).
- `pr_table`: pressure-ratio lookup table over (`n_corr_grid`, `w_corr_grid`).
- `eta_table`: isentropic-efficiency lookup table over (`n_corr_grid`, `w_corr_grid`).
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

"""
Corrected shaft speed from physical speed and local total temperature.
"""
corrected_speed(N::Real, Tt_in::Real, Tt_ref::Real) =
    N / sqrt(Tt_in / Tt_ref)

"""
Corrected mass flow from physical flow and local total conditions.
"""
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

"""
Interpolate pressure ratio and isentropic efficiency at corrected coordinates.
"""
function map_pr_eta(map::PerformanceMap, N_corr::Real, W_corr::Real)
    n = Float64(N_corr)
    w = Float64(W_corr)
    PR = _bilinear(n, w, map.n_corr_grid, map.w_corr_grid, map.pr_table)
    eta = _bilinear(n, w, map.n_corr_grid, map.w_corr_grid, map.eta_table)
    return (PR=PR, eta=eta)
end

"""
Compute corrected coordinates and return `(N_corr, W_corr, PR, eta)` from
physical shaft speed, mass flow, and local stagnation conditions.

Parameters:
- `map`: performance map wrapper.
- `N`: physical shaft speed.
- `mdot`: physical mass flow rate.
- `Tt_in`: local total (stagnation) temperature at the component inlet station.
- `Pt_in`: local total (stagnation) pressure at the component inlet station.
"""
function map_pr_eta_from_stagnation(
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

"""
Demo compressor map for testing and examples.

This is synthetic data intended for development/demo use only.
"""
function demo_compressor_map()
    PerformanceMap(
        288.15,
        101_325.0,
        [0.6, 0.8, 1.0],
        [12.0, 16.0, 20.0],
        [
            1.35 1.55 1.70;
            1.55 1.80 2.00;
            1.70 2.00 2.25;
        ],
        [
            0.74 0.78 0.75;
            0.78 0.83 0.80;
            0.76 0.82 0.79;
        ],
    )
end

"""
Demo turbine map for testing and examples.

This is synthetic data intended for development/demo use only.
Pressure-ratio convention here is `Pt_out / Pt_in` (less than 1 for expansion).
"""
function demo_turbine_map()
    PerformanceMap(
        288.15,
        101_325.0,
        [0.6, 0.8, 1.0],
        [10.0, 14.0, 18.0],
        [
            0.78 0.72 0.68;
            0.74 0.68 0.64;
            0.70 0.65 0.60;
        ],
        [
            0.82 0.84 0.83;
            0.85 0.88 0.87;
            0.84 0.87 0.86;
        ],
    )
end
