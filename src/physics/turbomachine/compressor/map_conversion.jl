"""
Compressor map conversion helpers between dimensional and non-dimensional tabulated forms.
"""

@inline function _linspace_inclusive(x_lo::Real, x_hi::Real, n::Int)
    n >= 2 || error("grid length must be at least 2")
    return collect(range(Float64(x_lo), Float64(x_hi), length=n))
end

@inline function _interp_linear_1d_clamped_conversion(
    xgrid::AbstractVector{<:Real},
    ygrid::AbstractVector{<:Real},
    x::Real,
)
    x1 = first(xgrid)
    x2 = last(xgrid)
    xc = clamp(x, x1, x2)
    i_hi = searchsortedfirst(xgrid, xc)
    i_hi <= 1 && return ygrid[1]
    i_hi > length(xgrid) && return ygrid[end]
    i_lo = i_hi - 1
    xl = xgrid[i_lo]
    xr = xgrid[i_hi]
    t = (xc - xl) / (xr - xl)
    return (1 - t) * ygrid[i_lo] + t * ygrid[i_hi]
end

"""
Convert a dimensional tabulated compressor map into a non-dimensional tabulated map.

This conversion is performed by resampling the source map onto a target
`(M_tip, phi_in)` grid at a chosen reference inlet state.

Arguments:
- `Tt_in_ref`, `Pt_in_ref`: inlet total reference state used for conversion.
- Geometry/gas: `gamma`, `gas_constant`, `tip_radius_inlet`, `mean_radius_inlet`, `inlet_area`.
- Optional target grids:
  - `m_tip_grid`: default is mapped from source `omega_corr_grid`.
  - `phi_in_grid`: default is linear over mapped source flow bounds.
"""
function to_nondimensional_tabulated_compressor_map(
    map::TabulatedCompressorPerformanceMap;
    gamma::Real,
    gas_constant::Real,
    tip_radius_inlet::Real,
    mean_radius_inlet::Real,
    inlet_area::Real,
    Tt_in_ref::Real=map.Tt_ref,
    Pt_in_ref::Real=map.Pt_ref,
    m_tip_grid::Union{Nothing,Vector{<:Real}}=nothing,
    phi_in_grid::Union{Nothing,Vector{<:Real}}=nothing,
    interpolation::Symbol=_interpolation_kind(map),
)
    Tt_in_ref > 0 || error("Tt_in_ref must be > 0")
    Pt_in_ref > 0 || error("Pt_in_ref must be > 0")
    gamma > 1 || error("gamma must be > 1")
    gas_constant > 0 || error("gas_constant must be > 0")
    tip_radius_inlet > 0 || error("tip_radius_inlet must be > 0")
    mean_radius_inlet > 0 || error("mean_radius_inlet must be > 0")
    inlet_area > 0 || error("inlet_area must be > 0")

    omega_grid_src = _omega_corr_grid(map)
    map_flow_grid_src = _mdot_corr_grid(map)
    surge_grid_src = map.mdot_corr_surge
    choke_grid_src = map.mdot_corr_choke
    all(omega_grid_src .> 0) || error("source omega grid must be strictly positive")

    a0_ref = sqrt(gamma * gas_constant * Tt_in_ref)
    m_grid = isnothing(m_tip_grid) ? (Float64.(omega_grid_src) .* Float64(tip_radius_inlet) ./ Float64(a0_ref)) : Float64.(m_tip_grid)
    length(m_grid) >= 2 || error("m_tip_grid must have at least 2 points")
    issorted(m_grid) || error("m_tip_grid must be sorted ascending")
    all(m_grid .> 0) || error("m_tip_grid values must be strictly positive")

    omega_from_m(m_tip) = Float64(m_tip) * Float64(a0_ref) / Float64(tip_radius_inlet)
    rho0_ref = Float64(Pt_in_ref) / (Float64(gas_constant) * Float64(Tt_in_ref))

    function _phi_from_source_map_flow(omega::Float64, map_flow::Float64)
        mdot = _physical_mdot_from_map_flow_coordinate(map, omega, map_flow, Tt_in_ref, Pt_in_ref)
        U_m1 = omega * Float64(mean_radius_inlet)
        return mdot / (rho0_ref * Float64(inlet_area) * U_m1)
    end

    phi_grid = if isnothing(phi_in_grid)
        phi_lo = minimum(
            _phi_from_source_map_flow(
                Float64(omega),
                _interp_linear_1d_clamped_conversion(omega_grid_src, surge_grid_src, Float64(omega)),
            ) for omega in omega_grid_src
        )
        phi_hi = maximum(
            _phi_from_source_map_flow(
                Float64(omega),
                _interp_linear_1d_clamped_conversion(omega_grid_src, choke_grid_src, Float64(omega)),
            ) for omega in omega_grid_src
        )
        _linspace_inclusive(phi_lo, phi_hi, length(map_flow_grid_src))
    else
        Float64.(phi_in_grid)
    end
    length(phi_grid) >= 2 || error("phi_in_grid must have at least 2 points")
    issorted(phi_grid) || error("phi_in_grid must be sorted ascending")

    pr = Matrix{Float64}(undef, length(m_grid), length(phi_grid))
    eta = Matrix{Float64}(undef, length(m_grid), length(phi_grid))
    for (i, m_tip) in pairs(m_grid)
        omega = omega_from_m(m_tip)
        for (j, phi) in pairs(phi_grid)
            mdot = Float64(phi) * rho0_ref * Float64(inlet_area) * omega * Float64(mean_radius_inlet)
            vals = performance_from_stagnation(map, omega, mdot, Tt_in_ref, Pt_in_ref)
            pr[i, j] = vals.PR
            eta[i, j] = vals.eta
        end
    end

    phi_surge = Float64[]
    phi_choke = Float64[]
    for m_tip in m_grid
        omega = omega_from_m(m_tip)
        map_flow_surge = _interp_linear_1d_clamped_conversion(omega_grid_src, surge_grid_src, omega)
        map_flow_choke = _interp_linear_1d_clamped_conversion(omega_grid_src, choke_grid_src, omega)
        push!(phi_surge, _phi_from_source_map_flow(omega, map_flow_surge))
        push!(phi_choke, _phi_from_source_map_flow(omega, map_flow_choke))
    end

    return NondimensionalPerformanceMap(
        gamma,
        gas_constant,
        tip_radius_inlet,
        mean_radius_inlet,
        inlet_area,
        m_grid,
        phi_grid,
        pr,
        eta;
        interpolation=interpolation,
        phi_surge=phi_surge,
        phi_choke=phi_choke,
    )
end

"""
Convert a non-dimensional tabulated compressor map into a dimensional tabulated map.

This conversion is performed by resampling the source map onto a target
`(omega, mdot_corr)` grid at a chosen reference inlet/reference-correction state.

Arguments:
- `Tt_in_ref`, `Pt_in_ref`: inlet total state used for physical conversion.
- `Tt_ref`, `Pt_ref`: corrected-flow reference of the output tabulated map.
- Optional target grids:
  - `omega_corr_grid`: default is mapped from source `m_tip_grid`.
  - `mdot_corr_grid`: default is linear over mapped source flow bounds.
"""
function to_tabulated_compressor_map(
    map::NondimensionalPerformanceMap;
    Tt_in_ref::Real=288.15,
    Pt_in_ref::Real=101_325.0,
    Tt_ref::Real=Tt_in_ref,
    Pt_ref::Real=Pt_in_ref,
    omega_corr_grid::Union{Nothing,Vector{<:Real}}=nothing,
    mdot_corr_grid::Union{Nothing,Vector{<:Real}}=nothing,
    interpolation::Symbol=_interpolation_kind(map),
)
    Tt_in_ref > 0 || error("Tt_in_ref must be > 0")
    Pt_in_ref > 0 || error("Pt_in_ref must be > 0")
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")

    m_grid_src = _m_tip_grid(map)
    phi_grid_src = _phi_in_grid(map)
    all(m_grid_src .> 0) || error("source M_tip grid must be strictly positive")

    a0_ref = sqrt(map.gamma * map.gas_constant * Tt_in_ref)
    omega_grid = isnothing(omega_corr_grid) ? (Float64.(m_grid_src) .* Float64(a0_ref) ./ map.tip_radius_inlet) : Float64.(omega_corr_grid)
    length(omega_grid) >= 2 || error("omega_corr_grid must have at least 2 points")
    issorted(omega_grid) || error("omega_corr_grid must be sorted ascending")
    all(omega_grid .> 0) || error("omega_corr_grid values must be strictly positive")

    rho0_ref = Float64(Pt_in_ref) / (map.gas_constant * Float64(Tt_in_ref))
    corr_fac = sqrt(Float64(Tt_in_ref) / Float64(Tt_ref)) / (Float64(Pt_in_ref) / Float64(Pt_ref))

    function _map_flow_from_source_phi(omega::Float64, phi::Float64)
        mdot = phi * rho0_ref * map.inlet_area * omega * map.mean_radius_inlet
        return mdot * corr_fac
    end

    mdot_grid = if isnothing(mdot_corr_grid)
        phi_lo = minimum(map.phi_surge)
        phi_hi = maximum(map.phi_choke)
        mdot_corr_lo = minimum(_map_flow_from_source_phi(Float64(omega), Float64(phi_lo)) for omega in omega_grid)
        mdot_corr_hi = maximum(_map_flow_from_source_phi(Float64(omega), Float64(phi_hi)) for omega in omega_grid)
        _linspace_inclusive(mdot_corr_lo, mdot_corr_hi, length(phi_grid_src))
    else
        Float64.(mdot_corr_grid)
    end
    length(mdot_grid) >= 2 || error("mdot_corr_grid must have at least 2 points")
    issorted(mdot_grid) || error("mdot_corr_grid must be sorted ascending")

    pr = Matrix{Float64}(undef, length(omega_grid), length(mdot_grid))
    eta = Matrix{Float64}(undef, length(omega_grid), length(mdot_grid))
    for (i, omega) in pairs(omega_grid)
        for (j, map_flow) in pairs(mdot_grid)
            mdot = Float64(map_flow) / corr_fac
            vals = performance_from_stagnation(map, omega, mdot, Tt_in_ref, Pt_in_ref)
            pr[i, j] = vals.PR
            eta[i, j] = vals.eta
        end
    end

    mdot_corr_surge = Float64[]
    mdot_corr_choke = Float64[]
    for omega in omega_grid
        m_tip = omega * map.tip_radius_inlet / a0_ref
        phi_s = _interp_linear_1d_clamped_conversion(m_grid_src, map.phi_surge, m_tip)
        phi_c = _interp_linear_1d_clamped_conversion(m_grid_src, map.phi_choke, m_tip)
        push!(mdot_corr_surge, _map_flow_from_source_phi(Float64(omega), phi_s))
        push!(mdot_corr_choke, _map_flow_from_source_phi(Float64(omega), phi_c))
    end

    return TabulatedCompressorPerformanceMap(
        Tt_ref,
        Pt_ref,
        omega_grid,
        mdot_grid,
        pr,
        eta;
        interpolation=interpolation,
        mdot_corr_surge=mdot_corr_surge,
        mdot_corr_choke=mdot_corr_choke,
    )
end
