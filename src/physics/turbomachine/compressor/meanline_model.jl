"""
Compressor meanline model wrappers built on shared axial-machine AxialMachine.
"""

import ..AxialMachine

const CompressorMeanlineModel = AxialMachine.AxialMachineModel

function _resolve_tabulation_grids(
    model::CompressorMeanlineModel,
    n_speed::Int,
    n_flow::Int,
    m_tip_grid::Union{Nothing,Vector{<:Real}},
    phi_in_grid::Union{Nothing,Vector{<:Real}},
)
    m_lo, m_hi = model.m_tip_bounds
    phi_lo, phi_hi = model.phi_in_bounds
    m_grid = isnothing(m_tip_grid) ? collect(range(m_lo, m_hi, length=n_speed)) : Float64.(m_tip_grid)
    phi_grid = isnothing(phi_in_grid) ? collect(range(phi_lo, phi_hi, length=n_flow)) : Float64.(phi_in_grid)
    return (m_grid=m_grid, phi_grid=phi_grid, phi_lo=phi_lo, phi_hi=phi_hi)
end

function _validate_tabulation_inputs(
    n_speed::Int,
    n_flow::Int,
    boundary_resolution::Int,
    interpolation::Symbol,
    m_grid::Vector{Float64},
    phi_grid::Vector{Float64},
)
    n_speed >= 2 || error("n_speed must be >= 2")
    n_flow >= 2 || error("n_flow must be >= 2")
    boundary_resolution >= 21 || error("boundary_resolution must be >= 21")
    interpolation in (:bilinear, :bicubic) ||
        error("interpolation must be :bilinear or :bicubic")
    length(m_grid) >= 2 || error("m_tip_grid must have at least 2 points")
    length(phi_grid) >= 2 || error("phi_in_grid must have at least 2 points")
    issorted(m_grid) || error("m_tip_grid must be sorted ascending")
    issorted(phi_grid) || error("phi_in_grid must be sorted ascending")
    return nothing
end

function _build_nd_tabulated_map(
    model::CompressorMeanlineModel,
    m_grid_valid::Vector{Float64},
    phi_grid::Vector{Float64},
    pr_table::Matrix{Float64},
    eta_table::Matrix{Float64},
    interpolation::Symbol,
    phi_surge::Vector{Float64},
    phi_choke::Vector{Float64},
)
    idx_ref = model.first_rotor_index
    tip_radius_inlet = model.rows[idx_ref].r_tip
    mean_radius_inlet = model.rows[idx_ref].r_mean
    inlet_area = model.A_ref * model.A_station[1]
    return NondimensionalPerformanceMap(
        model.gamma,
        model.gas_constant,
        tip_radius_inlet,
        mean_radius_inlet,
        inlet_area,
        m_grid_valid,
        phi_grid,
        pr_table,
        eta_table;
        interpolation=interpolation,
        phi_surge=phi_surge,
        phi_choke=phi_choke,
    )
end

"""
Tabulate a meanline compressor model onto a non-dimensional performance-map grid.
"""
function tabulate_compressor_meanline_model(
    model::CompressorMeanlineModel;
    n_speed::Int=31,
    n_flow::Int=41,
    m_tip_grid::Union{Nothing,Vector{<:Real}}=nothing,
    phi_in_grid::Union{Nothing,Vector{<:Real}}=nothing,
    interpolation::Symbol=:bilinear,
    boundary_resolution::Int=401,
)
    grids = _resolve_tabulation_grids(model, n_speed, n_flow, m_tip_grid, phi_in_grid)
    _validate_tabulation_inputs(
        n_speed,
        n_flow,
        boundary_resolution,
        interpolation,
        grids.m_grid,
        grids.phi_grid,
    )

    limits = AxialMachine.feasible_flow_limits(
        model,
        grids.m_grid,
        grids.phi_lo,
        grids.phi_hi;
        boundary_resolution=boundary_resolution,
        is_feasible=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta)),
    )
    length(limits.valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid speed lines in requested tabulation range")
    m_grid_valid = grids.m_grid[limits.valid_speed_idx]
    phi_surge = limits.flow_min
    phi_choke = limits.flow_max
    tables = AxialMachine.sample_streamtube_solve(
        model,
        m_grid_valid,
        grids.phi_grid;
        flow_min=phi_surge,
        flow_max=phi_choke,
        is_feasible=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta)),
    )
    return _build_nd_tabulated_map(
        model,
        m_grid_valid,
        grids.phi_grid,
        tables.pr_table,
        tables.eta_table,
        interpolation,
        phi_surge,
        phi_choke,
    )
end

"""
Demo compressor meanline model for development/testing.
"""
demo_compressor_meanline_model() = AxialMachine.demo_axial_machine_model()
