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
    tip_radius_inlet = model.r_tip_ref
    mean_radius_inlet = abs(model.speed_ratio_ref) * model.r_flow_ref
    inlet_area = AxialMachine.station_area(model, 1)
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

function _resolve_dim_tabulation_grids(
    model::CompressorMeanlineModel,
    n_speed::Int,
    n_flow::Int,
    omega_corr_grid::Union{Nothing,Vector{<:Real}},
    mdot_corr_grid::Union{Nothing,Vector{<:Real}},
    Tt_in_ref::Real,
    Pt_in_ref::Real,
    Tt_ref::Real,
    Pt_ref::Real,
    omega_ref_for_phi::Union{Nothing,Real},
)
    n_speed >= 2 || error("n_speed must be >= 2")
    n_flow >= 2 || error("n_flow must be >= 2")
    Tt_in_ref > 0 || error("Tt_in_ref must be > 0")
    Pt_in_ref > 0 || error("Pt_in_ref must be > 0")
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")

    a0_in_ref = sqrt(model.gamma * model.gas_constant * Tt_in_ref)
    m_lo, m_hi = model.m_tip_bounds
    omega_lo = m_lo * a0_in_ref / model.r_tip_ref
    omega_hi = m_hi * a0_in_ref / model.r_tip_ref
    omega_grid = isnothing(omega_corr_grid) ?
        collect(range(omega_lo, omega_hi, length=n_speed)) :
        Float64.(omega_corr_grid)
    length(omega_grid) >= 2 || error("omega_corr_grid must have at least 2 points")
    issorted(omega_grid) || error("omega_corr_grid must be sorted ascending")
    all(omega_grid .> 0) || error("omega_corr_grid values must be strictly positive")

    omega_ref_for_phi_f = isnothing(omega_ref_for_phi) ?
        0.5 * (first(omega_grid) + last(omega_grid)) :
        Float64(omega_ref_for_phi)
    omega_ref_for_phi_f > 0 || error("omega_ref_for_phi must be > 0")

    mean_radius_inlet = abs(model.speed_ratio_ref) * model.r_flow_ref
    inlet_area = AxialMachine.station_area(model, 1)
    rho0_in_ref = Pt_in_ref / (model.gas_constant * Tt_in_ref)
    corr_fac = sqrt(Tt_in_ref / Tt_ref) / (Pt_in_ref / Pt_ref)
    phi_lo, phi_hi = model.phi_in_bounds
    mdot_corr_lo = phi_lo * rho0_in_ref * inlet_area * omega_ref_for_phi_f * mean_radius_inlet * corr_fac
    mdot_corr_hi = phi_hi * rho0_in_ref * inlet_area * omega_ref_for_phi_f * mean_radius_inlet * corr_fac

    mdot_corr_grid_f = isnothing(mdot_corr_grid) ?
        collect(range(mdot_corr_lo, mdot_corr_hi, length=n_flow)) :
        Float64.(mdot_corr_grid)
    length(mdot_corr_grid_f) >= 2 || error("mdot_corr_grid must have at least 2 points")
    issorted(mdot_corr_grid_f) || error("mdot_corr_grid must be sorted ascending")

    return (
        omega_corr_grid=omega_grid,
        mdot_corr_grid=mdot_corr_grid_f,
        omega_ref_for_phi=omega_ref_for_phi_f,
    )
end

"""
Tabulate a meanline compressor model onto a non-dimensional performance-map grid.
"""
function tabulate_compressor_meanline_model_nd(
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
Tabulate a meanline compressor model onto a dimensional corrected-flow performance-map grid.

The output map is built on a uniform `(omega_corr, mdot_corr)` grid (unless explicit
grids are provided), while streamtube sampling is executed in non-dimensional coordinates.
"""
function tabulate_compressor_meanline_model_dim(
    model::CompressorMeanlineModel;
    n_speed::Int=31,
    n_flow::Int=41,
    omega_corr_grid::Union{Nothing,Vector{<:Real}}=nothing,
    mdot_corr_grid::Union{Nothing,Vector{<:Real}}=nothing,
    Tt_in_ref::Real=288.15,
    Pt_in_ref::Real=101_325.0,
    Tt_ref::Real=Tt_in_ref,
    Pt_ref::Real=Pt_in_ref,
    omega_ref_for_phi::Union{Nothing,Real}=nothing,
    interpolation::Symbol=:bilinear,
    boundary_resolution::Int=401,
)
    boundary_resolution >= 21 || error("boundary_resolution must be >= 21")
    interpolation in (:bilinear, :bicubic) ||
        error("interpolation must be :bilinear or :bicubic")

    dim_grids = _resolve_dim_tabulation_grids(
        model,
        n_speed,
        n_flow,
        omega_corr_grid,
        mdot_corr_grid,
        Tt_in_ref,
        Pt_in_ref,
        Tt_ref,
        Pt_ref,
        omega_ref_for_phi,
    )

    physical_grids = corrected_grids_to_physical_grids(
        dim_grids.omega_corr_grid,
        dim_grids.mdot_corr_grid;
        Tt_in=Tt_in_ref,
        Pt_in=Pt_in_ref,
        Tt_ref=Tt_ref,
        Pt_ref=Pt_ref,
    )

    mean_radius_inlet = abs(model.speed_ratio_ref) * model.r_flow_ref
    inlet_area = AxialMachine.station_area(model, 1)
    nd_grids = physical_grids_to_nondimensional_grids(
        physical_grids.omega_grid,
        physical_grids.mdot_grid;
        gamma=model.gamma,
        gas_constant=model.gas_constant,
        tip_radius_inlet=model.r_tip_ref,
        mean_radius_inlet=mean_radius_inlet,
        inlet_area=inlet_area,
        Tt_in=Tt_in_ref,
        Pt_in=Pt_in_ref,
        omega_ref_for_phi=dim_grids.omega_ref_for_phi,
    )

    limits = AxialMachine.feasible_flow_limits(
        model,
        nd_grids.m_tip_grid,
        minimum(nd_grids.phi_in_grid),
        maximum(nd_grids.phi_in_grid);
        boundary_resolution=boundary_resolution,
        is_feasible=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta)),
    )
    length(limits.valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid speed lines in requested tabulation range")

    m_grid_valid = nd_grids.m_tip_grid[limits.valid_speed_idx]
    omega_corr_grid_valid = dim_grids.omega_corr_grid[limits.valid_speed_idx]
    phi_surge = limits.flow_min
    phi_choke = limits.flow_max

    tables = AxialMachine.sample_streamtube_solve(
        model,
        m_grid_valid,
        nd_grids.phi_in_grid;
        flow_min=phi_surge,
        flow_max=phi_choke,
        is_feasible=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta)),
    )

    nd_map = _build_nd_tabulated_map(
        model,
        m_grid_valid,
        nd_grids.phi_in_grid,
        tables.pr_table,
        tables.eta_table,
        interpolation,
        phi_surge,
        phi_choke,
    )
    return to_tabulated_compressor_map(
        nd_map;
        Tt_in_ref=Tt_in_ref,
        Pt_in_ref=Pt_in_ref,
        Tt_ref=Tt_ref,
        Pt_ref=Pt_ref,
        omega_corr_grid=omega_corr_grid_valid,
        mdot_corr_grid=dim_grids.mdot_corr_grid,
        interpolation=interpolation,
    )
end

"""
Demo compressor meanline model for development/testing.
"""
demo_compressor_meanline_model() = AxialMachine.demo_axial_compressor_model()

"""
Demo turbine-like meanline model serialized via compressor meanline schema.
"""
demo_turbine_meanline_model() = AxialMachine.demo_axial_turbine_model()
