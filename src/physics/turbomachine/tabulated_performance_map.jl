"""
Tabulated performance map implementation.
"""

abstract type AbstractPerformanceMap end

"""
Tabulated compressor/turbine performance map on corrected coordinates.

Fields:
- `Tt_ref`: reference total temperature for corrected normalization.
- `Pt_ref`: reference total pressure for corrected normalization.
- `omega_corr_grid`: corrected speed grid (ascending).
- `mdot_corr_grid`: corrected mass-flow grid (ascending).
- `pr_table`: total-pressure-ratio table (`Pt_out/Pt_in`).
- `eta_table`: adiabatic-efficiency table.
"""
struct TabulatedPerformanceMap <: AbstractPerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    omega_corr_grid::Vector{Float64}
    mdot_corr_grid::Vector{Float64}
    pr_table::Matrix{Float64}
    eta_table::Matrix{Float64}
end

function TabulatedPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    omega_corr_grid::Vector{<:Real},
    mdot_corr_grid::Vector{<:Real},
    pr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real},
)
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(omega_corr_grid) >= 2 || error("omega_corr_grid must have at least 2 points")
    length(mdot_corr_grid) >= 2 || error("mdot_corr_grid must have at least 2 points")
    issorted(omega_corr_grid) || error("omega_corr_grid must be sorted ascending")
    issorted(mdot_corr_grid) || error("mdot_corr_grid must be sorted ascending")
    size(pr_table) == (length(omega_corr_grid), length(mdot_corr_grid)) ||
        error("pr_table size must match (length(omega_corr_grid), length(mdot_corr_grid))")
    size(eta_table) == (length(omega_corr_grid), length(mdot_corr_grid)) ||
        error("eta_table size must match (length(omega_corr_grid), length(mdot_corr_grid))")

    return TabulatedPerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        Float64.(omega_corr_grid),
        Float64.(mdot_corr_grid),
        Float64.(pr_table),
        Float64.(eta_table),
    )
end

"""Corrected shaft speed from physical speed and local total temperature."""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) =
    omega / sqrt(Tt_in / Tt_ref)

"""Corrected mass flow from physical flow and local total conditions."""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(omega::Real, Tt_in::Real, map::AbstractPerformanceMap) =
    corrected_speed(omega, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::AbstractPerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

@inline _map_primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x

function _cell_index_and_fraction(grid::AbstractVector{<:Real}, x::Real)
    x_primal = _map_primal_value(x)
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

function _bilinear(
    x::Real,
    y::Real,
    xgrid::AbstractVector{<:Real},
    ygrid::AbstractVector{<:Real},
    table::AbstractMatrix{<:Real},
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
Evaluate a turbomachine map at corrected coordinates.

Returns named tuple `(PR, eta)` where:
- `PR` is total-pressure ratio (`Pt_out/Pt_in`)
- `eta` is adiabatic efficiency
"""
function performance_map(
    map::TabulatedPerformanceMap,
    omega_corr::Real,
    mdot_corr::Real,
)
    PR = _bilinear(
        omega_corr,
        mdot_corr,
        map.omega_corr_grid,
        map.mdot_corr_grid,
        map.pr_table,
    )
    eta = _bilinear(
        omega_corr,
        mdot_corr,
        map.omega_corr_grid,
        map.mdot_corr_grid,
        map.eta_table,
    )
    return (PR=PR, eta=eta)
end

"""
Evaluate a turbomachine map from physical values and local stagnation state.

Returns `(omega_corr, mdot_corr, PR, eta)`.
"""
function performance_map_from_stagnation(
    map::AbstractPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    omega_corr = corrected_speed(omega, Tt_in, map)
    mdot_corr = corrected_flow(mdot, Tt_in, Pt_in, map)
    vals = performance_map(map, omega_corr, mdot_corr)
    return (omega_corr=omega_corr, mdot_corr=mdot_corr, PR=vals.PR, eta=vals.eta)
end

"""Demo tabulated compressor map for development/testing."""
function demo_compressor_performance_map()
    TabulatedPerformanceMap(
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
Demo tabulated turbine map for development/testing.
Pressure-ratio convention is `Pt_out / Pt_in` (< 1 for expansion).
"""
function demo_turbine_performance_map()
    TabulatedPerformanceMap(
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
