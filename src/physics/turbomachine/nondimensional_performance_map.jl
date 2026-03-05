using ...Utility: AbstractTableMap, interpolation_map, table_xgrid, table_ygrid

"""
Shared tabulated non-dimensional turbomachine performance map.

Domain inputs:
- `M_tip = U_tip / a0_in`
- `phi_in = Vx_1 / U_m_1`

Outputs:
- `PR = Pt_out / Pt_in`
- `eta` total-to-total efficiency

Geometry fields are used to convert physical `(omega, mdot, Tt_in, Pt_in)` into
`(M_tip, phi_in)`.
"""
struct NondimensionalPerformanceMap{M<:AbstractTableMap} <: AbstractCompressorPerformanceMap
    gamma::Float64
    gas_constant::Float64
    tip_radius_inlet::Float64
    mean_radius_inlet::Float64
    inlet_area::Float64
    pr_map::M
    eta_map::M
    phi_surge::Vector{Float64}
    phi_choke::Vector{Float64}
end

function NondimensionalPerformanceMap(
    gamma::Real,
    gas_constant::Real,
    tip_radius_inlet::Real,
    mean_radius_inlet::Real,
    inlet_area::Real,
    pr_map::M,
    eta_map::M,
    phi_surge::Vector{<:Real},
    phi_choke::Vector{<:Real},
) where {M<:AbstractTableMap}
    gamma > 1 || error("gamma must be > 1")
    gas_constant > 0 || error("gas_constant must be > 0")
    tip_radius_inlet > 0 || error("tip_radius_inlet must be > 0")
    mean_radius_inlet > 0 || error("mean_radius_inlet must be > 0")
    inlet_area > 0 || error("inlet_area must be > 0")
    table_xgrid(pr_map) == table_xgrid(eta_map) || error("pr_map/eta_map x grids must match")
    table_ygrid(pr_map) == table_ygrid(eta_map) || error("pr_map/eta_map y grids must match")

    m_grid = table_xgrid(pr_map)
    length(phi_surge) == length(m_grid) || error("phi_surge length must match M_tip grid length")
    length(phi_choke) == length(m_grid) || error("phi_choke length must match M_tip grid length")

    phi_surge_f = Float64.(phi_surge)
    phi_choke_f = Float64.(phi_choke)
    all(phi_surge_f .<= phi_choke_f) || error("phi_surge must be <= phi_choke at every M_tip grid point")

    return NondimensionalPerformanceMap(
        Float64(gamma),
        Float64(gas_constant),
        Float64(tip_radius_inlet),
        Float64(mean_radius_inlet),
        Float64(inlet_area),
        pr_map,
        eta_map,
        phi_surge_f,
        phi_choke_f,
    )
end

function NondimensionalPerformanceMap(
    gamma::Real,
    gas_constant::Real,
    tip_radius_inlet::Real,
    mean_radius_inlet::Real,
    inlet_area::Real,
    m_tip_grid::Vector{<:Real},
    phi_in_grid::Vector{<:Real},
    pr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real};
    interpolation::Symbol,
    phi_surge::Union{Nothing,Vector{<:Real}}=nothing,
    phi_choke::Union{Nothing,Vector{<:Real}}=nothing,
)
    length(m_tip_grid) >= 2 || error("m_tip_grid must have at least 2 points")
    length(phi_in_grid) >= 2 || error("phi_in_grid must have at least 2 points")
    issorted(m_tip_grid) || error("m_tip_grid must be sorted ascending")
    issorted(phi_in_grid) || error("phi_in_grid must be sorted ascending")
    size(pr_table) == (length(m_tip_grid), length(phi_in_grid)) ||
        error("pr_table size must match (length(m_tip_grid), length(phi_in_grid))")
    size(eta_table) == (length(m_tip_grid), length(phi_in_grid)) ||
        error("eta_table size must match (length(m_tip_grid), length(phi_in_grid))")

    m_tip_grid_f = Float64.(m_tip_grid)
    phi_in_grid_f = Float64.(phi_in_grid)
    pr_table_f = Float64.(pr_table)
    eta_table_f = Float64.(eta_table)
    pr_map = interpolation_map(interpolation, m_tip_grid_f, phi_in_grid_f, pr_table_f)
    eta_map = interpolation_map(interpolation, m_tip_grid_f, phi_in_grid_f, eta_table_f)

    phi_surge_f = isnothing(phi_surge) ? fill(first(phi_in_grid_f), length(m_tip_grid_f)) : Float64.(phi_surge)
    phi_choke_f = isnothing(phi_choke) ? fill(last(phi_in_grid_f), length(m_tip_grid_f)) : Float64.(phi_choke)

    return NondimensionalPerformanceMap(
        gamma,
        gas_constant,
        tip_radius_inlet,
        mean_radius_inlet,
        inlet_area,
        pr_map,
        eta_map,
        phi_surge_f,
        phi_choke_f,
    )
end
