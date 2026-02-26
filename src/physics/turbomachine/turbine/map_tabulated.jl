"""
Tabulated turbine performance map implementation.
"""

using ....Utility: BilinearMap, bilinear_evaluate

abstract type AbstractTurbinePerformanceMap end

"""
Tabulated turbine performance map on corrected coordinates.

Fields:
- `Tt_ref`: reference total temperature for corrected normalization.
- `Pt_ref`: reference total pressure for corrected normalization.
- `omega_corr_grid`: corrected speed grid (ascending).
- `pr_turb_grid`: turbine expansion-ratio grid (ascending), where
  `PR_turb = Pt_in / Pt_out`.
- `mdot_corr_table`: corrected mass-flow table.
- `eta_table`: adiabatic-efficiency table.
"""
struct TabulatedTurbinePerformanceMap <: AbstractTurbinePerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    omega_corr_grid::Vector{Float64}
    pr_turb_grid::Vector{Float64}
    mdot_corr_table::Matrix{Float64}
    eta_table::Matrix{Float64}
    mdot_corr_map::BilinearMap
    eta_map::BilinearMap
end

function TabulatedTurbinePerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    omega_corr_grid::Vector{<:Real},
    pr_turb_grid::Vector{<:Real},
    mdot_corr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real},
)
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(omega_corr_grid) >= 2 || error("omega_corr_grid must have at least 2 points")
    length(pr_turb_grid) >= 2 || error("pr_turb_grid must have at least 2 points")
    issorted(omega_corr_grid) || error("omega_corr_grid must be sorted ascending")
    issorted(pr_turb_grid) || error("pr_turb_grid must be sorted ascending")
    size(mdot_corr_table) == (length(omega_corr_grid), length(pr_turb_grid)) ||
        error("mdot_corr_table size must match (length(omega_corr_grid), length(pr_turb_grid))")
    size(eta_table) == (length(omega_corr_grid), length(pr_turb_grid)) ||
        error("eta_table size must match (length(omega_corr_grid), length(pr_turb_grid))")

    omega_corr_grid_f = Float64.(omega_corr_grid)
    pr_turb_grid_f = Float64.(pr_turb_grid)
    mdot_corr_table_f = Float64.(mdot_corr_table)
    eta_table_f = Float64.(eta_table)
    mdot_corr_map = BilinearMap(omega_corr_grid_f, pr_turb_grid_f, mdot_corr_table_f)
    eta_map = BilinearMap(omega_corr_grid_f, pr_turb_grid_f, eta_table_f)

    return TabulatedTurbinePerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        omega_corr_grid_f,
        pr_turb_grid_f,
        mdot_corr_table_f,
        eta_table_f,
        mdot_corr_map,
        eta_map,
    )
end

"""Corrected shaft speed from physical speed and local total temperature."""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) =
    omega / sqrt(Tt_in / Tt_ref)

"""Corrected mass flow from physical flow and local total conditions."""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(omega::Real, Tt_in::Real, map::AbstractTurbinePerformanceMap) =
    corrected_speed(omega, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::AbstractTurbinePerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

_physical_flow_from_corrected(
    mdot_corr::Real,
    Tt_in::Real,
    Pt_in::Real,
    Tt_ref::Real,
    Pt_ref::Real,
) = mdot_corr * (Pt_in / Pt_ref) / sqrt(Tt_in / Tt_ref)

"""
Evaluate a turbine map at corrected coordinates.

Returns named tuple `(mdot_corr, eta)` where:
- `mdot_corr` is corrected mass flow
- `eta` is adiabatic efficiency
"""
function turbine_performance_map(
    map::TabulatedTurbinePerformanceMap,
    omega_corr::Real,
    pr_turb::Real,
)
    mdot_corr = bilinear_evaluate(map.mdot_corr_map, omega_corr, pr_turb)
    eta = bilinear_evaluate(map.eta_map, omega_corr, pr_turb)
    return (mdot_corr=mdot_corr, eta=eta)
end

"""
Evaluate a turbine map from physical values and local stagnation state.

Returns `(omega_corr, PR_turb, mdot_corr, mdot, eta)` where
`PR_turb = Pt_in / Pt_out`.
"""
function turbine_performance_map_from_stagnation(
    map::AbstractTurbinePerformanceMap,
    omega::Real,
    Pt_in::Real,
    Pt_out::Real,
    Tt_in::Real,
)
    Pt_out > 0 || error("Pt_out must be > 0")
    PR_turb = Pt_in / Pt_out
    omega_corr = corrected_speed(omega, Tt_in, map)
    vals = turbine_performance_map(map, omega_corr, PR_turb)
    mdot = _physical_flow_from_corrected(vals.mdot_corr, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)
    return (
        omega_corr=omega_corr,
        PR_turb=PR_turb,
        mdot_corr=vals.mdot_corr,
        mdot=mdot,
        eta=vals.eta,
    )
end

"""Demo tabulated turbine map for development/testing."""
function demo_turbine_performance_map()
    TabulatedTurbinePerformanceMap(
        288.15,
        101_325.0,
        [0.6, 0.8, 1.0],
        [1.4, 1.8, 2.2],
        [
            10.0 12.0 13.5;
            12.5 14.5 16.0;
            14.0 16.5 18.5;
        ],
        [
            0.80 0.84 0.82;
            0.83 0.87 0.85;
            0.81 0.86 0.84;
        ],
    )
end
