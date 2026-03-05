"""
Compressor meanline model wrappers built on shared axial-machine core.
"""

import ..AxialMachine

const Core = AxialMachine
const CompressorMeanlineModel = Core.AxialMachineModel

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

function _compute_surge_choke_limits(
    model::CompressorMeanlineModel,
    m_grid::Vector{Float64},
    phi_lo::Float64,
    phi_hi::Float64,
    boundary_resolution::Int,
)
    phi_probe = collect(range(phi_lo, phi_hi, length=boundary_resolution))
    valid_speed_idx = Int[]
    phi_surge = Float64[]
    phi_choke = Float64[]
    for (i, m_tip) in pairs(m_grid)
        valid_phis = Float64[]
        for phi in phi_probe
            vals = Core.streamtube_solve(model, m_tip, phi)
            if vals.valid && isfinite(vals.PR) && isfinite(vals.eta)
                push!(valid_phis, phi)
            end
        end
        isempty(valid_phis) && continue
        push!(valid_speed_idx, i)
        push!(phi_surge, first(valid_phis))
        push!(phi_choke, last(valid_phis))
    end
    length(valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid speed lines in requested tabulation range")
    return (
        m_grid_valid=m_grid[valid_speed_idx],
        phi_surge=phi_surge,
        phi_choke=phi_choke,
    )
end

function _compute_nd_tables(
    model::CompressorMeanlineModel,
    m_grid_valid::Vector{Float64},
    phi_grid::Vector{Float64},
    phi_surge::Vector{Float64},
    phi_choke::Vector{Float64},
)
    pr_table = Matrix{Float64}(undef, length(m_grid_valid), length(phi_grid))
    eta_table = Matrix{Float64}(undef, length(m_grid_valid), length(phi_grid))
    for (i, m_tip) in pairs(m_grid_valid)
        for (j, phi_raw) in pairs(phi_grid)
            phi = clamp(phi_raw, phi_surge[i], phi_choke[i])
            vals = Core.streamtube_solve(model, m_tip, phi)
            (vals.valid && isfinite(vals.PR) && isfinite(vals.eta)) ||
                error("meanline tabulation produced non-finite value at m_tip=$(m_tip), phi=$(phi)")
            pr_table[i, j] = vals.PR
            eta_table[i, j] = vals.eta
        end
    end
    return (pr_table=pr_table, eta_table=eta_table)
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
    idx_ref = Core.first_rotor_index(model)
    tip_radius_inlet = model.rows[idx_ref].r_tip
    mean_radius_inlet = model.rows[idx_ref].r_mean
    inlet_area = model.A_ref * model.A_station[1]
    return NonDimensionalTabulatedCompressorPerformanceMap(
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

    limits = _compute_surge_choke_limits(
        model,
        grids.m_grid,
        grids.phi_lo,
        grids.phi_hi,
        boundary_resolution,
    )
    tables = _compute_nd_tables(
        model,
        limits.m_grid_valid,
        grids.phi_grid,
        limits.phi_surge,
        limits.phi_choke,
    )
    return _build_nd_tabulated_map(
        model,
        limits.m_grid_valid,
        grids.phi_grid,
        tables.pr_table,
        tables.eta_table,
        interpolation,
        limits.phi_surge,
        limits.phi_choke,
    )
end

"""
Demo compressor meanline model for development/testing.
"""
demo_compressor_meanline_model() = Core.demo_axial_machine_model()
