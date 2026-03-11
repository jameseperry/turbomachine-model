"""
Shared axial-machine meanline tabulation.

This file builds the canonical `TabulatedPerformanceMap` directly from the
axial-machine streamtube solver. The resulting map stores corrected
`(omega_corr, mdot_corr)` tables internally, while callers continue to evaluate
it using physical `(omega, mdot)` through `performance_from_stagnation(...)`.
"""

using ..Fluids

const AxialModel = AxialMachine.AxialMachineModel
const AxialMeanlineModel = AxialModel

function _default_axial_eos(model::AxialModel)
    return Fluids.IdealGasEOS(:axial; gas_constant=model.gas_constant, gamma=model.gamma)
end

demo_axial_compressor_model() = AxialMachine.demo_axial_compressor_model()
demo_axial_turbine_model() = AxialMachine.demo_axial_turbine_model()

function _resolve_meanline_tabulation_grids(
    model::AxialModel,
    eos::Fluids.AbstractEOS,
    n_speed::Int,
    n_flow::Int,
    omega_corr_grid::Union{Nothing,AbstractVector{<:Real}},
    mdot_corr_grid::Union{Nothing,AbstractVector{<:Real}},
    Tt_in_ref::Real,
    Pt_in_ref::Real,
    Tt_ref::Real,
    Pt_ref::Real,
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

    ht_in_ref = Fluids.enthalpy_from_temperature(eos, Float64(Tt_in_ref))
    a0_in_ref = Fluids.speed_of_sound(eos, Float64(Pt_in_ref), ht_in_ref)
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

    inlet_area = AxialMachine.station_area(model, 1)
    Vtheta_inlet_ref = 0.0
    radii_f = Float64.(streamtube_radii)
    _mdot_from_Vx(Vx) = AxialMachine._build_station_state(
        eos,
        AxialMachine._station_radius(model, radii_f, 1),
        inlet_area,
        Float64(Pt_in_ref),
        ht_in_ref,
        Vx,
        Vtheta_inlet_ref,
        NaN,
        true,
        false,
        false,
    ).mdot_station

    if isnothing(mdot_corr_grid)
        Vx_grid_f = collect(range(model.Vx_bounds[1], model.Vx_bounds[2], length=n_flow))
        mdot_grid = [_mdot_from_Vx(Vx) for Vx in Vx_grid_f]
        all(isfinite, mdot_grid) || error("Vx_bounds produce invalid inlet mass-flow values at the reference state")
        issorted(mdot_grid) || error("projected mdot grid from Vx_bounds must be sorted ascending")
        mdot_corr_grid_f = corrected_flow.(
            mdot_grid,
            Float64(Tt_in_ref),
            Float64(Pt_in_ref),
            Float64(Tt_ref),
            Float64(Pt_ref),
        )
    else
        mdot_corr_grid_f = Float64.(mdot_corr_grid)
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
        Vx_grid_f = Float64[]
        for mdot in mdot_grid
            inlet = Fluids.velocity_from_stagnation_massflow(
                eos,
                mdot,
                Float64(Pt_in_ref),
                ht_in_ref,
                inlet_area,
                Vtheta_inlet_ref;
                prefer=:low,
            )
            inlet.converged || error("mdot_corr_grid contains an inlet mass flow that does not admit a valid Vx solution")
            push!(Vx_grid_f, inlet.Vx)
        end
    end

    length(mdot_corr_grid_f) >= 2 || error("mdot_corr_grid must have at least 2 points")
    issorted(mdot_corr_grid_f) || error("mdot_corr_grid must be sorted ascending")
    all(mdot_corr_grid_f .> 0) || error("mdot_corr_grid values must be strictly positive")

    return (
        omega_corr_grid=omega_corr_grid_f,
        omega_grid=omega_grid,
        Vx_grid=Vx_grid_f,
        mdot_corr_grid=mdot_corr_grid_f,
        mdot_grid=mdot_grid,
        ht_in_ref=ht_in_ref,
        a0_in_ref=a0_in_ref,
        streamtube_radii=radii_f,
    )
end

function _sample_axial_machine_streamtube(
    model::AxialModel,
    eos::Fluids.AbstractEOS,
    omega::Real,
    Vx_inlet::Real;
    pt_in::Real,
    ht_in::Real,
    streamtube_radii::AbstractVector{<:Real},
    Vtheta_inlet::Real,
    prefer_root::Symbol,
)
    vals = AxialMachine.streamtube_solve(
        model,
        eos,
        streamtube_radii,
        Float64(pt_in),
        Float64(ht_in),
        Float64(omega),
        Float64(Vx_inlet),
        Float64(Vtheta_inlet);
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
    eos::Fluids.AbstractEOS,
    grids::NamedTuple,
    Tt_in_ref::Real,
    Pt_in_ref::Real,
    Tt_ref::Real,
    Pt_ref::Real;
    pt_in::Real,
    ht_in::Real,
    boundary_resolution::Int,
    Vtheta_inlet::Real,
    prefer_root::Symbol,
)
    Vx_probe = collect(range(model.Vx_bounds[1], model.Vx_bounds[2], length=boundary_resolution))
    valid_speed_idx = Int[]
    Vx_min = Float64[]
    Vx_max = Float64[]
    mdot_min = Float64[]
    mdot_max = Float64[]

    for (i, omega) in pairs(grids.omega_grid)
        feasible_Vx = Float64[]
        feasible_mdot = Float64[]
        for Vx in Vx_probe
            vals = _sample_axial_machine_streamtube(
                model,
                eos,
                omega,
                Vx;
                pt_in=Float64(pt_in),
                ht_in=Float64(ht_in),
                streamtube_radii=grids.streamtube_radii,
                Vtheta_inlet=Float64(Vtheta_inlet),
                prefer_root=prefer_root,
            )
            if vals.valid && isfinite(vals.PR) && isfinite(vals.eta) && vals.PR > 0 && isfinite(vals.raw.stations[1].mdot_station)
                push!(feasible_Vx, Vx)
                push!(feasible_mdot, vals.raw.stations[1].mdot_station)
            end
        end
        isempty(feasible_Vx) && continue
        push!(valid_speed_idx, i)
        push!(Vx_min, first(feasible_Vx))
        push!(Vx_max, last(feasible_Vx))
        push!(mdot_min, first(feasible_mdot))
        push!(mdot_max, last(feasible_mdot))
    end

    length(valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid speed lines in requested tabulation range")

    mdot_corr_min = Float64[]
    mdot_corr_max = Float64[]
    for i in eachindex(valid_speed_idx)
        push!(mdot_corr_min, corrected_flow(mdot_min[i], Float64(Tt_in_ref), Float64(Pt_in_ref), Float64(Tt_ref), Float64(Pt_ref)))
        push!(mdot_corr_max, corrected_flow(mdot_max[i], Float64(Tt_in_ref), Float64(Pt_in_ref), Float64(Tt_ref), Float64(Pt_ref)))
    end

    return (
        valid_speed_idx=valid_speed_idx,
        Vx_min=Vx_min,
        Vx_max=Vx_max,
        mdot_min=mdot_min,
        mdot_max=mdot_max,
        mdot_corr_min=mdot_corr_min,
        mdot_corr_max=mdot_corr_max,
    )
end

function _sample_meanline_tables(
    model::AxialModel,
    eos::Fluids.AbstractEOS,
    grids::NamedTuple,
    limits::NamedTuple;
    pt_in::Real,
    ht_in::Real,
    Vtheta_inlet::Real,
    prefer_root::Symbol,
    want_diagnostics::Bool,
)
    valid_idx = limits.valid_speed_idx
    omega_corr_valid = grids.omega_corr_grid[valid_idx]
    omega_valid = grids.omega_grid[valid_idx]
    n_speed = length(valid_idx)
    n_flow = length(grids.mdot_corr_grid)
    pr_table = Matrix{Float64}(undef, n_speed, n_flow)
    eta_table = Matrix{Float64}(undef, n_speed, n_flow)
    diagnostics = want_diagnostics ? Matrix{NamedTuple}(undef, n_speed, n_flow) : Matrix{NamedTuple}(undef, 0, 0)

    for i in eachindex(valid_idx)
        for j in 1:n_flow
            Vx_used = clamp(grids.Vx_grid[j], limits.Vx_min[i], limits.Vx_max[i])
            sample = _sample_axial_machine_streamtube(
                model,
                eos,
                omega_valid[i],
                Vx_used;
                pt_in=Float64(pt_in),
                ht_in=Float64(ht_in),
                streamtube_radii=grids.streamtube_radii,
                Vtheta_inlet=Float64(Vtheta_inlet),
                prefer_root=prefer_root,
            )
            pr_table[i, j] = sample.PR
            eta_table[i, j] = sample.eta

            if want_diagnostics
                diagnostics[i, j] = (
                    omega_corr=omega_corr_valid[i],
                    omega=omega_valid[i],
                    Vx=grids.Vx_grid[j],
                    Vx_used=Vx_used,
                    mdot_corr=grids.mdot_corr_grid[j],
                    mdot=grids.mdot_grid[j],
                    mdot_used=sample.raw.stations[1].mdot_station,
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
    eos::Fluids.AbstractEOS=_default_axial_eos(model),
    n_speed::Int=31,
    n_flow::Int=41,
    omega_corr_grid::Union{Nothing,AbstractVector{<:Real}}=nothing,
    mdot_corr_grid::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Tt_in_ref::Real=288.15,
    Pt_in_ref::Real=101_325.0,
    Tt_ref::Real=Tt_in_ref,
    Pt_ref::Real=Pt_in_ref,
    interpolation::Symbol=:bilinear,
    boundary_resolution::Int=401,
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
    want_diagnostics::Bool=true,
)
    interpolation in (:bilinear, :bicubic) ||
        error("interpolation must be :bilinear or :bicubic")
    boundary_resolution >= 21 || error("boundary_resolution must be >= 21")

    grids = _resolve_meanline_tabulation_grids(
        model,
        eos,
        n_speed,
        n_flow,
        omega_corr_grid,
        mdot_corr_grid,
        Tt_in_ref,
        Pt_in_ref,
        Tt_ref,
        Pt_ref,
        streamtube_radii,
    )
    limits = _compute_feasible_flow_limits(
        model,
        eos,
        grids,
        Tt_in_ref,
        Pt_in_ref,
        Tt_ref,
        Pt_ref;
        pt_in=Float64(Pt_in_ref),
        ht_in=grids.ht_in_ref,
        boundary_resolution=boundary_resolution,
        Vtheta_inlet=Float64(Vtheta_inlet),
        prefer_root=prefer_root,
    )
    tables = _sample_meanline_tables(
        model,
        eos,
        grids,
        limits;
        pt_in=Float64(Pt_in_ref),
        ht_in=grids.ht_in_ref,
        Vtheta_inlet=Float64(Vtheta_inlet),
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
