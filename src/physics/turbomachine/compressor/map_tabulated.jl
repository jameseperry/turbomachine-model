"""
Tabulated compressor performance map implementation.
"""

using ....Utility: BilinearMap, bilinear_evaluate

abstract type AbstractCompressorPerformanceMap end

"""
Tabulated compressor performance map on corrected coordinates.

Fields:
- `Tt_ref`: reference total temperature for corrected normalization.
- `Pt_ref`: reference total pressure for corrected normalization.
- `omega_corr_grid`: corrected speed grid (ascending).
- `mdot_corr_grid`: corrected mass-flow grid (ascending).
- `pr_table`: total-pressure-ratio table (`Pt_out/Pt_in`).
- `eta_table`: adiabatic-efficiency table.
"""
struct TabulatedCompressorPerformanceMap <: AbstractCompressorPerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    omega_corr_grid::Vector{Float64}
    mdot_corr_grid::Vector{Float64}
    pr_table::Matrix{Float64}
    eta_table::Matrix{Float64}
    pr_map::BilinearMap
    eta_map::BilinearMap
end

function TabulatedCompressorPerformanceMap(
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

    omega_corr_grid_f = Float64.(omega_corr_grid)
    mdot_corr_grid_f = Float64.(mdot_corr_grid)
    pr_table_f = Float64.(pr_table)
    eta_table_f = Float64.(eta_table)
    pr_map = BilinearMap(omega_corr_grid_f, mdot_corr_grid_f, pr_table_f)
    eta_map = BilinearMap(omega_corr_grid_f, mdot_corr_grid_f, eta_table_f)

    return TabulatedCompressorPerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        omega_corr_grid_f,
        mdot_corr_grid_f,
        pr_table_f,
        eta_table_f,
        pr_map,
        eta_map,
    )
end

"""Corrected shaft speed from physical speed and local total temperature."""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) =
    omega / sqrt(Tt_in / Tt_ref)

"""Corrected mass flow from physical flow and local total conditions."""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(omega::Real, Tt_in::Real, map::AbstractCompressorPerformanceMap) =
    corrected_speed(omega, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::AbstractCompressorPerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

"""
Evaluate a compressor map at corrected coordinates.

Returns named tuple `(PR, eta)` where:
- `PR` is total-pressure ratio (`Pt_out/Pt_in`)
- `eta` is adiabatic efficiency
"""
function compressor_performance_map(
    map::TabulatedCompressorPerformanceMap,
    omega_corr::Real,
    mdot_corr::Real,
)
    PR = bilinear_evaluate(map.pr_map, omega_corr, mdot_corr)
    eta = bilinear_evaluate(map.eta_map, omega_corr, mdot_corr)
    return (PR=PR, eta=eta)
end

"""
Evaluate a compressor map from physical values and local stagnation state.

Returns `(omega_corr, mdot_corr, PR, eta)`.
"""
function compressor_performance_map_from_stagnation(
    map::AbstractCompressorPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    omega_corr = corrected_speed(omega, Tt_in, map)
    mdot_corr = corrected_flow(mdot, Tt_in, Pt_in, map)
    vals = compressor_performance_map(map, omega_corr, mdot_corr)
    return (omega_corr=omega_corr, mdot_corr=mdot_corr, PR=vals.PR, eta=vals.eta)
end

"""Demo tabulated compressor map for development/testing."""
function demo_compressor_performance_map()
    TabulatedCompressorPerformanceMap(
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
