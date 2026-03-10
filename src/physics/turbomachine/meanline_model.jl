"""
Shared axial-machine meanline tabulation.

This file builds the canonical `TabulatedPerformanceMap` directly from the
axial-machine streamtube solver. The resulting map stores corrected
`(omega_corr, mdot_corr)` tables internally, while callers continue to evaluate
it using physical `(omega, mdot)` through `performance_from_stagnation(...)`.
"""

const AxialModel = AxialMachine.AxialMachineModel
const AxialMeanlineModel = AxialModel

demo_axial_compressor_model() = AxialMachine.demo_axial_compressor_model()
demo_axial_turbine_model() = AxialMachine.demo_axial_turbine_model()

function _resolve_meanline_tabulation_grids(
    model::AxialModel,
    n_speed::Int,
    n_flow::Int,
    omega_corr_grid::Union{Nothing,AbstractVector{<:Real}},
    mdot_corr_grid::Union{Nothing,AbstractVector{<:Real}},
    Tt_in_ref::Real,
    Pt_in_ref::Real,
    Tt_ref::Real,
    Pt_ref::Real,
    omega_ref_for_phi::Union{Nothing,Real},
    streamtube_radii::AbstractVector{<:Real},
)
    n_speed >= 2 || error("n_speed must be >= 2")
    n_flow >= 2 || error("n_flow must be >= 2")
    Tt_in_ref > 0 || error("Tt_in_ref must be > 0")
    Pt_in_ref > 0 || error("Pt_in_ref must be > 0")
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(streamtube_radii) == length(model.rows) ||
        error("streamtube_radii length must match number of rows")

    a0_in_ref = sqrt(model.gamma * model.gas_constant * Float64(Tt_in_ref))
    m_lo, m_hi = model.m_tip_bounds
    omega_lo = m_lo * a0_in_ref / model.r_tip_ref
    omega_hi = m_hi * a0_in_ref / model.r_tip_ref
    omega_corr_lo = corrected_speed(omega_lo, Tt_in_ref, Tt_ref)
    omega_corr_hi = corrected_speed(omega_hi, Tt_in_ref, Tt_ref)
    omega_corr_grid_f = isnothing(omega_corr_grid) ?
        collect(range(omega_corr_lo, omega_corr_hi, length=n_speed)) :
        Float64.(omega_corr_grid)
    length(omega_corr_grid_f) >= 2 || error("omega_corr_grid must have at least 2 points")
    issorted(omega_corr_grid_f) || error("omega_corr_grid must be sorted ascending")
    all(omega_corr_grid_f .> 0) || error("omega_corr_grid values must be strictly positive")

    omega_grid = physical_speed.(omega_corr_grid_f, Float64(Tt_in_ref), Float64(Tt_ref))
    omega_ref_for_phi_f = isnothing(omega_ref_for_phi) ?
        0.5 * (first(omega_grid) + last(omega_grid)) :
        Float64(omega_ref_for_phi)
    omega_ref_for_phi_f > 0 || error("omega_ref_for_phi must be > 0")

    mean_radius_inlet = abs(model.speed_ratio_ref) * model.r_flow_ref
    inlet_area = AxialMachine.station_area(model, 1)
    rho0_in_ref = Float64(Pt_in_ref) / (model.gas_constant * Float64(Tt_in_ref))
    corr_fac = sqrt(Float64(Tt_in_ref) / Float64(Tt_ref)) / (Float64(Pt_in_ref) / Float64(Pt_ref))
    phi_lo, phi_hi = model.phi_in_bounds
    mdot_corr_lo = phi_lo * rho0_in_ref * inlet_area * omega_ref_for_phi_f * mean_radius_inlet * corr_fac
    mdot_corr_hi = phi_hi * rho0_in_ref * inlet_area * omega_ref_for_phi_f * mean_radius_inlet * corr_fac
    mdot_corr_grid_f = isnothing(mdot_corr_grid) ?
        collect(range(mdot_corr_lo, mdot_corr_hi, length=n_flow)) :
        Float64.(mdot_corr_grid)
    length(mdot_corr_grid_f) >= 2 || error("mdot_corr_grid must have at least 2 points")
    issorted(mdot_corr_grid_f) || error("mdot_corr_grid must be sorted ascending")
    all(mdot_corr_grid_f .> 0) || error("mdot_corr_grid values must be strictly positive")

    mdot_grid = physical_flow.(
        mdot_corr_grid_f,
        Float64(Tt_in_ref),
        Float64(Pt_in_ref),
        Float64(Tt_ref),
        Float64(Pt_ref),
    )
    m_tip_grid = omega_grid .* model.r_tip_ref ./ a0_in_ref

    return (
        omega_corr_grid=omega_corr_grid_f,
        omega_grid=omega_grid,
        mdot_corr_grid=mdot_corr_grid_f,
        mdot_grid=mdot_grid,
        m_tip_grid=m_tip_grid,
        rho0_in_ref=rho0_in_ref,
        inlet_area=inlet_area,
        mean_radius_inlet=mean_radius_inlet,
        corr_fac=corr_fac,
        a0_in_ref=a0_in_ref,
        omega_ref_for_phi=omega_ref_for_phi_f,
        streamtube_radii=Float64.(streamtube_radii),
    )
end

function _sample_axial_machine_streamtube(
    model::AxialModel,
    m_tip::Real,
    phi_in::Real;
    streamtube_radii::AbstractVector{<:Real},
    nu_theta_inlet::Real,
    prefer_root::Symbol,
)
    vals = AxialMachine.streamtube_solve_with_phi(
        model,
        m_tip,
        phi_in;
        streamtube_radii=streamtube_radii,
        nu_theta_inlet=nu_theta_inlet,
        prefer_root=prefer_root,
    )
    return (
        PR=Float64(vals.PR),
        eta=Float64(vals.eta),
        valid=Bool(vals.valid),
        stall=Bool(vals.stall),
        choke=Bool(vals.choke),
        raw=vals,
    )
end

function _compute_feasible_flow_limits(
    model::AxialModel,
    grids::NamedTuple;
    boundary_resolution::Int,
    nu_theta_inlet::Real,
    prefer_root::Symbol,
)
    limits = AxialMachine.feasible_flow_limits(
        model,
        grids.m_tip_grid,
        model.phi_in_bounds[1],
        model.phi_in_bounds[2];
        boundary_resolution=boundary_resolution,
        streamtube_radii=grids.streamtube_radii,
        nu_theta_inlet=nu_theta_inlet,
        prefer_root=prefer_root,
        is_feasible=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta) && vals.PR > 0),
    )

    length(limits.valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid speed lines in requested tabulation range")

    mdot_corr_min = Float64[]
    mdot_corr_max = Float64[]
    for idx in limits.valid_speed_idx
        omega = grids.omega_grid[idx]
        scale = grids.rho0_in_ref * grids.inlet_area * abs(omega * grids.mean_radius_inlet) * grids.corr_fac
        push!(mdot_corr_min, limits.flow_min[length(mdot_corr_min) + 1] * scale)
        push!(mdot_corr_max, limits.flow_max[length(mdot_corr_max) + 1] * scale)
    end

    return (
        valid_speed_idx=limits.valid_speed_idx,
        phi_min=limits.flow_min,
        phi_max=limits.flow_max,
        mdot_corr_min=mdot_corr_min,
        mdot_corr_max=mdot_corr_max,
    )
end

function _sample_meanline_tables(
    model::AxialModel,
    grids::NamedTuple,
    limits::NamedTuple;
    nu_theta_inlet::Real,
    prefer_root::Symbol,
    want_diagnostics::Bool,
)
    valid_idx = limits.valid_speed_idx
    omega_corr_valid = grids.omega_corr_grid[valid_idx]
    omega_valid = grids.omega_grid[valid_idx]
    m_tip_valid = grids.m_tip_grid[valid_idx]
    n_speed = length(valid_idx)
    n_flow = length(grids.mdot_corr_grid)
    pr_table = Matrix{Float64}(undef, n_speed, n_flow)
    eta_table = Matrix{Float64}(undef, n_speed, n_flow)
    diagnostics = want_diagnostics ? Matrix{NamedTuple}(undef, n_speed, n_flow) : Matrix{NamedTuple}(undef, 0, 0)

    for i in eachindex(valid_idx)
        omega_i = omega_valid[i]
        denom = grids.rho0_in_ref * grids.inlet_area * abs(omega_i * grids.mean_radius_inlet)
        phi_targets = grids.mdot_grid ./ denom
        row_tables = AxialMachine.sample_streamtube_solve(
            model,
            [m_tip_valid[i]],
            phi_targets;
            flow_min=[limits.phi_min[i]],
            flow_max=[limits.phi_max[i]],
            streamtube_radii=grids.streamtube_radii,
            nu_theta_inlet=nu_theta_inlet,
            prefer_root=prefer_root,
            is_feasible=(vals -> vals.valid && isfinite(vals.PR) && isfinite(vals.eta) && vals.PR > 0),
        )
        pr_table[i, :] = row_tables.pr_table[1, :]
        eta_table[i, :] = row_tables.eta_table[1, :]

        if want_diagnostics
            for j in 1:n_flow
                phi_target = phi_targets[j]
                phi_used = clamp(phi_target, limits.phi_min[i], limits.phi_max[i])
                sample = _sample_axial_machine_streamtube(
                    model,
                    m_tip_valid[i],
                    phi_used;
                    streamtube_radii=grids.streamtube_radii,
                    nu_theta_inlet=nu_theta_inlet,
                    prefer_root=prefer_root,
                )
                diagnostics[i, j] = (
                    omega_corr=omega_corr_valid[i],
                    omega=omega_i,
                    mdot_corr=grids.mdot_corr_grid[j],
                    mdot=grids.mdot_grid[j],
                    phi_target=phi_target,
                    phi_used=phi_used,
                    m_tip=m_tip_valid[i],
                    PR=sample.PR,
                    eta=sample.eta,
                    valid=sample.valid,
                    stall=sample.stall,
                    choke=sample.choke,
                )
            end
        end
    end

    return (
        omega_corr_valid=omega_corr_valid,
        pr_table=pr_table,
        eta_table=eta_table,
        diagnostics=diagnostics,
    )
end

function _build_tabulated_performance_map(
    grids::NamedTuple,
    limits::NamedTuple,
    tables::NamedTuple,
    Tt_ref::Real,
    Pt_ref::Real,
    interpolation::Symbol,
)
    return TabulatedPerformanceMap(
        Tt_ref,
        Pt_ref,
        tables.omega_corr_valid,
        grids.mdot_corr_grid,
        tables.pr_table,
        tables.eta_table;
        interpolation=interpolation,
        flow_min=limits.mdot_corr_min,
        flow_max=limits.mdot_corr_max,
    )
end

"""
    tabulate_axial_machine_model(model; kwargs...)

Tabulate a shared `AxialMachineModel` into the canonical turbomachine
`TabulatedPerformanceMap`.

The output map stores corrected `(omega_corr, mdot_corr)` tables internally,
while runtime callers continue to use physical `(omega, mdot)` through the
shared map API.
"""
function tabulate_axial_machine_model(
    model::AxialModel;
    n_speed::Int=31,
    n_flow::Int=41,
    omega_corr_grid::Union{Nothing,AbstractVector{<:Real}}=nothing,
    mdot_corr_grid::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Tt_in_ref::Real=288.15,
    Pt_in_ref::Real=101_325.0,
    Tt_ref::Real=Tt_in_ref,
    Pt_ref::Real=Pt_in_ref,
    omega_ref_for_phi::Union{Nothing,Real}=nothing,
    interpolation::Symbol=:bilinear,
    boundary_resolution::Int=401,
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    nu_theta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
    want_diagnostics::Bool=true,
)
    interpolation in (:bilinear, :bicubic) ||
        error("interpolation must be :bilinear or :bicubic")
    boundary_resolution >= 21 || error("boundary_resolution must be >= 21")

    grids = _resolve_meanline_tabulation_grids(
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
        streamtube_radii,
    )
    limits = _compute_feasible_flow_limits(
        model,
        grids;
        boundary_resolution=boundary_resolution,
        nu_theta_inlet=nu_theta_inlet,
        prefer_root=prefer_root,
    )
    tables = _sample_meanline_tables(
        model,
        grids,
        limits;
        nu_theta_inlet=nu_theta_inlet,
        prefer_root=prefer_root,
        want_diagnostics=want_diagnostics,
    )
    return _build_tabulated_performance_map(
        grids,
        limits,
        tables,
        Tt_ref,
        Pt_ref,
        interpolation,
    )
end

"""
Compatibility alias while callers migrate to `tabulate_axial_machine_model`.
"""
function tabulate_axial_machine_model_dim(model::AxialModel; kwargs...)
    return tabulate_axial_machine_model(model; kwargs...)
end

const CompressorMeanlineModel = AxialModel
const TurbineMeanlineModel = AxialModel

function tabulate_compressor_meanline_model_dim(model::AxialModel; kwargs...)
    return tabulate_axial_machine_model(model; kwargs...)
end

function tabulate_turbine_meanline_model(model::AxialModel; kwargs...)
    return tabulate_axial_machine_model(model; kwargs...)
end
