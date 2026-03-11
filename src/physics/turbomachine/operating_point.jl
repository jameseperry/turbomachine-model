"""
Generic operating-point solve for the canonical `TabulatedPerformanceMap`.

The solver is branch-complete: it returns every physically admissible `mdot`
branch that satisfies the requested `(pt_in, ht_in, pt_out, omega)` operating
condition, and leaves branch selection to downstream code.
"""

using ..Fluids
using ...Utility: bracket_bisect_roots

@inline function _safe_omega_divisor(omega::Float64)
    abs(omega) > 1e-12 || return copysign(1e-12, omega == 0.0 ? 1.0 : omega)
    return omega
end

@inline function _safe_enthalpy_from_temperature(
    eos::Fluids.AbstractEOS,
    T::Float64,
)
    return (isfinite(T) && T > 0.0) ? Fluids.enthalpy_from_temperature(eos, T) : NaN
end

function _operating_point_mdot_bounds(
    map::TabulatedPerformanceMap,
    omega::Float64,
    Tt_in::Float64,
    pt_in::Float64,
)
    domain = performance_map_domain(map, Tt_in, pt_in)
    mdot_lo = domain.mdot_flow_range.min(omega)
    mdot_hi = domain.mdot_flow_range.max(omega)
    if !(isfinite(mdot_lo) && isfinite(mdot_hi))
        return (NaN, NaN)
    end
    return (min(mdot_lo, mdot_hi), max(mdot_lo, mdot_hi))
end

function _reconstruct_ht_out(
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    pt_out::Float64,
    eta::Float64,
)
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    if h2s >= ht_in
        eta > 0.0 || return NaN
        return ht_in + (h2s - ht_in) / eta
    end
    eta > 0.0 || return NaN
    return ht_in - eta * (ht_in - h2s)
end

function _operating_point_residuals(
    pt_in::Float64,
    ht_in::Float64,
    pt_out::Float64,
    ht_out::Float64,
    mdot::Float64,
    omega::Float64,
    tau::Float64,
    PR::Float64,
    eta::Float64,
    h2s::Float64,
)
    h_is = h2s - ht_in
    h_actual = ht_out - ht_in
    R_pout = pt_out - PR * pt_in
    R_eff = if h_is >= 0.0
        eta * h_actual - h_is
    else
        h_actual - eta * h_is
    end
    R_P = tau * omega - mdot * h_actual
    return (R_pout, R_eff, R_P)
end

function _is_physically_admissible(
    eos::Fluids.AbstractEOS,
    pt_out::Float64,
    ht_out::Float64,
    mdot::Float64,
    PR::Float64,
    eta::Float64,
    tau::Float64,
    power::Float64,
)
    if !(isfinite(pt_out) && isfinite(ht_out) && isfinite(mdot) &&
         isfinite(PR) && isfinite(eta) && isfinite(tau) && isfinite(power))
        return false
    end
    pt_out > 0.0 || return false
    ht_out > 0.0 || return false
    try
        Tt_out = Fluids.temperature(eos, pt_out, ht_out)
        return isfinite(Tt_out) && Tt_out > 0.0
    catch
        return false
    end
end

function _thermo_efficiency(
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    pt_out::Float64,
    ht_out::Float64,
)
    if !(isfinite(pt_in) && isfinite(ht_in) && isfinite(pt_out) && isfinite(ht_out))
        return NaN
    end
    if !(pt_in > 0.0 && pt_out > 0.0)
        return NaN
    end
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    if h2s >= ht_in
        denom = ht_out - ht_in
        abs(denom) > 1e-12 || return NaN
        return (h2s - ht_in) / denom
    end
    denom = ht_in - h2s
    abs(denom) > 1e-12 || return NaN
    return (ht_in - ht_out) / denom
end

function _package_streamtube_diagnostics(
    raw::AxialMachine.StreamtubeSolveResult,
    inputs::NamedTuple,
    summary::NamedTuple,
)
    return (
        inputs=inputs,
        physical=(
            station_data=raw.station_data,
            row_data=raw.row_data,
        ),
        summary=summary,
    )
end

function _inlet_Vx_from_mdot(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    mdot::Float64,
    Vtheta_inlet::Float64,
    prefer_root::Symbol,
)
    inlet = Fluids.velocity_from_stagnation_massflow(
        eos,
        mdot,
        pt_in,
        ht_in,
        AxialMachine.station_area(model, 1),
        Vtheta_inlet;
        prefer=prefer_root,
    )
    return inlet.converged ? inlet.Vx : NaN
end

function _stage_diagnostics_from_candidate(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS,
    candidate::NamedTuple;
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    streamtube_radii::AbstractVector{<:Real},
    Vtheta_inlet::Float64,
    prefer_root::Symbol,
)
    Vx_inlet = _inlet_Vx_from_mdot(
        model,
        eos,
        pt_in,
        ht_in,
        Float64(candidate.mdot),
        Vtheta_inlet,
        prefer_root,
    )
    raw = AxialMachine.streamtube_solve(
        model,
        eos,
        streamtube_radii,
        pt_in,
        ht_in,
        omega,
        Vx_inlet,
        Vtheta_inlet;
        prefer_root=prefer_root,
    )
    ht_out_streamtube = raw.ht_t[end]
    tau_streamtube = Float64(candidate.mdot) * (ht_out_streamtube - ht_in) / _safe_omega_divisor(omega)
    power_streamtube = tau_streamtube * omega
    operating_point_thermo_efficiency = _thermo_efficiency(
        eos,
        pt_in,
        ht_in,
        Float64(candidate.pressure_ratio) * pt_in,
        Float64(candidate.ht_out),
    )
    streamtube_thermo_efficiency = _thermo_efficiency(
        eos,
        pt_in,
        ht_in,
        raw.pressure_ratio * pt_in,
        ht_out_streamtube,
    )
    return _package_streamtube_diagnostics(
        raw,
        (
            pt_in=pt_in,
            ht_in=ht_in,
            omega=omega,
            mdot=raw.stations[1].mdot_station,
            Vx_inlet=raw.stations[1].Vx,
            Vtheta_inlet=Vtheta_inlet,
            streamtube_radii=Float64.(streamtube_radii),
        ),
        (
            operating_point_pressure_ratio=Float64(candidate.pressure_ratio),
            operating_point_efficiency=Float64(candidate.efficiency),
            operating_point_thermo_efficiency=operating_point_thermo_efficiency,
            operating_point_ht_out=Float64(candidate.ht_out),
            operating_point_tau=Float64(candidate.tau),
            operating_point_power=Float64(candidate.power),
            streamtube_pressure_ratio=raw.pressure_ratio,
            streamtube_efficiency=raw.efficiency,
            streamtube_thermo_efficiency=streamtube_thermo_efficiency,
            streamtube_ht_out=ht_out_streamtube,
            streamtube_tau=tau_streamtube,
            streamtube_power=power_streamtube,
            pressure_ratio_error=raw.pressure_ratio - Float64(candidate.pressure_ratio),
            efficiency_error=raw.efficiency - Float64(candidate.efficiency),
            thermo_efficiency_error=streamtube_thermo_efficiency - operating_point_thermo_efficiency,
            ht_out_error=ht_out_streamtube - Float64(candidate.ht_out),
            tau_error=tau_streamtube - Float64(candidate.tau),
            power_error=power_streamtube - Float64(candidate.power),
            valid=raw.valid,
            stall=raw.stall,
            choke=raw.choke,
        ),
    )
end

"""
    diagnose_axial_operating_point(model, eos; pt_in, ht_in, omega, Vx_inlet, streamtube_radii=meanline_radii(model), Vtheta_inlet=0.0, prefer_root=:low)

Run the axial-machine streamtube solver directly at one physical operating point
and return detailed row/station data from the dimensional solver.
"""
function diagnose_axial_operating_point(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    Vx_inlet::Real,
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
)
    pt_in_f = Float64(pt_in)
    ht_in_f = Float64(ht_in)
    omega_f = Float64(omega)
    Vx_inlet_f = Float64(Vx_inlet)
    raw = AxialMachine.streamtube_solve(
        model,
        eos,
        streamtube_radii,
        pt_in_f,
        ht_in_f,
        omega_f,
        Vx_inlet_f,
        Float64(Vtheta_inlet);
        prefer_root=prefer_root,
    )
    Tt_out = raw.station_data[end].Tt
    ht_out = raw.ht_t[end]
    mdot_f = raw.stations[1].mdot_station
    tau_shaft = mdot_f * (ht_out - ht_in_f) / _safe_omega_divisor(omega_f)
    power = tau_shaft * omega_f
    thermo_efficiency = _thermo_efficiency(eos, pt_in_f, ht_in_f, raw.pressure_ratio * pt_in_f, ht_out)
    return _package_streamtube_diagnostics(
        raw,
        (
            pt_in=pt_in_f,
            ht_in=ht_in_f,
            omega=omega_f,
            mdot=mdot_f,
            Vx_inlet=raw.stations[1].Vx,
            Vtheta_inlet=Float64(Vtheta_inlet),
            streamtube_radii=Float64.(streamtube_radii),
        ),
        (
            mdot=mdot_f,
            omega=omega_f,
            pressure_ratio=raw.pressure_ratio,
            efficiency=raw.efficiency,
            thermo_efficiency=thermo_efficiency,
            Tt_out=Tt_out,
            ht_out=ht_out,
            tau=tau_shaft,
            power=power,
            valid=raw.valid,
            stall=raw.stall,
            choke=raw.choke,
        ),
    )
end

"""
    replay_operating_point_with_streamtube(model, eos, candidate; pt_in, ht_in, omega, streamtube_radii=meanline_radii(model), Vtheta_inlet=0.0, prefer_root=:low)

Re-run the source axial-machine streamtube model at a solved operating-point
branch and return detailed stage diagnostics.

This is a post-processing tool for comparing the tabulated-map operating point
to the causal row-marching model that generated the map. The returned summary
includes both the map branch values and the replayed streamtube values, along
with the mismatch between them.
"""
function replay_operating_point_with_streamtube(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS,
    candidate::NamedTuple;
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
)
    return _stage_diagnostics_from_candidate(
        model,
        eos,
        candidate;
        pt_in=Float64(pt_in),
        ht_in=Float64(ht_in),
        omega=Float64(omega),
        streamtube_radii=Float64.(streamtube_radii),
        Vtheta_inlet=Float64(Vtheta_inlet),
        prefer_root=prefer_root,
    )
end

function replay_operating_point_with_streamtube(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    Vx_inlet::Real,
    pressure_ratio::Real,
    efficiency::Real,
    ht_out::Real,
    tau::Real,
    power::Real=(Float64(tau) * Float64(omega)),
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    Vtheta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
)
    pt_in_f = Float64(pt_in)
    ht_in_f = Float64(ht_in)
    omega_f = Float64(omega)
    Vx_inlet_f = Float64(Vx_inlet)
    Vtheta_inlet_f = Float64(Vtheta_inlet)
    radii_f = Float64.(streamtube_radii)
    raw = AxialMachine.streamtube_solve(
        model,
        eos,
        radii_f,
        pt_in_f,
        ht_in_f,
        omega_f,
        Vx_inlet_f,
        Vtheta_inlet_f;
        prefer_root=prefer_root,
    )
    ht_out_streamtube = raw.ht_t[end]
    tau_streamtube = raw.stations[1].mdot_station * (ht_out_streamtube - ht_in_f) / _safe_omega_divisor(omega_f)
    power_streamtube = tau_streamtube * omega_f
    operating_point_thermo_efficiency = _thermo_efficiency(
        eos,
        pt_in_f,
        ht_in_f,
        Float64(pressure_ratio) * pt_in_f,
        Float64(ht_out),
    )
    streamtube_thermo_efficiency = _thermo_efficiency(
        eos,
        pt_in_f,
        ht_in_f,
        raw.pressure_ratio * pt_in_f,
        ht_out_streamtube,
    )
    return _package_streamtube_diagnostics(
        raw,
        (
            pt_in=pt_in_f,
            ht_in=ht_in_f,
            omega=omega_f,
            mdot=raw.stations[1].mdot_station,
            Vx_inlet=raw.stations[1].Vx,
            Vtheta_inlet=Vtheta_inlet_f,
            streamtube_radii=radii_f,
        ),
        (
            operating_point_pressure_ratio=Float64(pressure_ratio),
            operating_point_efficiency=Float64(efficiency),
            operating_point_thermo_efficiency=operating_point_thermo_efficiency,
            operating_point_ht_out=Float64(ht_out),
            operating_point_tau=Float64(tau),
            operating_point_power=Float64(power),
            streamtube_pressure_ratio=raw.pressure_ratio,
            streamtube_efficiency=raw.efficiency,
            streamtube_thermo_efficiency=streamtube_thermo_efficiency,
            streamtube_ht_out=ht_out_streamtube,
            streamtube_tau=tau_streamtube,
            streamtube_power=power_streamtube,
            pressure_ratio_error=raw.pressure_ratio - Float64(pressure_ratio),
            efficiency_error=raw.efficiency - Float64(efficiency),
            thermo_efficiency_error=streamtube_thermo_efficiency - operating_point_thermo_efficiency,
            ht_out_error=ht_out_streamtube - Float64(ht_out),
            tau_error=tau_streamtube - Float64(tau),
            power_error=power_streamtube - Float64(power),
            valid=raw.valid,
            stall=raw.stall,
            choke=raw.choke,
        ),
    )
end

function _operating_point_candidates(
    map::TabulatedPerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Float64,
    ht_in::Float64,
    pt_out::Float64,
    omega::Float64,
    continuation_hints::AbstractVector{<:Real}=Float64[],
    want_diagnostics::Bool=true,
    options::NamedTuple=NamedTuple(),
)
    pt_in > 0.0 || error("pt_in must be > 0")
    pt_out > 0.0 || error("pt_out must be > 0")

    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    target_pr = pt_out / pt_in
    mdot_lo, mdot_hi = _operating_point_mdot_bounds(map, omega, Tt_in, pt_in)
    n_scan = get(options, :n_scan, 401)
    pr_tol = Float64(get(options, :pr_tol, 1e-8))
    continuation_band_fraction = Float64(get(options, :continuation_band_fraction, 0.02))
    max_bisect_iters = Int(get(options, :max_bisect_iters, 60))

    if !(isfinite(mdot_lo) && isfinite(mdot_hi) && (mdot_hi > mdot_lo))
        diagnostics = (
            pt_in=pt_in,
            pt_out=pt_out,
            ht_in=ht_in,
            omega=omega,
            target_pr=target_pr,
            Tt_in=Tt_in,
            mdot_bounds=(mdot_lo, mdot_hi),
            reason=:degenerate_mdot_bounds,
        )
        return (candidates=NamedTuple[], rejected_candidates=NamedTuple[], diagnostics=diagnostics)
    end

    residual = mdot -> begin
        vals = performance_from_stagnation(map, omega, mdot, Tt_in, pt_in)
        return vals.pressure_ratio - target_pr
    end

    roots = bracket_bisect_roots(
        residual,
        (mdot_lo, mdot_hi);
        n_scan=n_scan,
        root_tol=pr_tol,
        prior_roots=continuation_hints,
        continuation_band_fraction=continuation_band_fraction,
        max_bisect_iters=max_bisect_iters,
    )

    candidates = NamedTuple[]
    rejected_candidates = NamedTuple[]
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    for mdot in roots
        vals = performance_from_stagnation(map, omega, mdot, Tt_in, pt_in)
        ht_out = _reconstruct_ht_out(eos, pt_in, ht_in, pt_out, vals.eta)
        tau = mdot * (ht_out - ht_in) / _safe_omega_divisor(omega)
        power = tau * omega
        residuals = _operating_point_residuals(
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau,
            vals.pressure_ratio,
            vals.eta,
            h2s,
        )
        candidate = (
            branch_coordinate=mdot,
            mdot=mdot,
            ht_out=ht_out,
            tau=tau,
            pressure_ratio=vals.pressure_ratio,
            efficiency=vals.eta,
            power=power,
            residuals=residuals,
            physically_admissible=_is_physically_admissible(
                eos,
                pt_out,
                ht_out,
                mdot,
                vals.pressure_ratio,
                vals.eta,
                tau,
                power,
            ),
            machine_payload=(
                speed_coord=vals.speed_coord,
                flow_coord=vals.flow_coord,
                low_flow=vals.low_flow,
                high_flow=vals.high_flow,
                valid=vals.valid,
            ),
        )
        if candidate.physically_admissible
            push!(candidates, candidate)
        elseif want_diagnostics
            push!(rejected_candidates, candidate)
        end
    end

    diagnostics = (
        omega=omega,
        pt_in=pt_in,
        pt_out=pt_out,
        ht_in=ht_in,
        Tt_in=Tt_in,
        target_pr=target_pr,
        mdot_bounds=(mdot_lo, mdot_hi),
        n_roots=length(roots),
        n_candidates=length(candidates),
        n_rejected_candidates=length(rejected_candidates),
        rejected_candidates=(want_diagnostics ? rejected_candidates : NamedTuple[]),
    )

    return (candidates=candidates, diagnostics=diagnostics)
end


"""
    solve_operating_point(map, eos; pt_in, ht_in, pt_out, omega, continuation_hints=Float64[], want_diagnostics=true, options=NamedTuple())

Solve a single operating point for the canonical `TabulatedPerformanceMap`.

Returns all physically admissible branches. Branch selection is intentionally
left to downstream code.
"""
function solve_operating_point(
    map::TabulatedPerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    omega::Real,
    continuation_hints::AbstractVector{<:Real}=Float64[],
    want_diagnostics::Bool=true,
    options::NamedTuple=NamedTuple(),
)
    point = _operating_point_candidates(
        map,
        eos;
        pt_in=Float64(pt_in),
        ht_in=Float64(ht_in),
        pt_out=Float64(pt_out),
        omega=Float64(omega),
        continuation_hints=continuation_hints,
        want_diagnostics=want_diagnostics,
        options=options,
    )
    return (
        converged=!isempty(point.candidates),
        retcode=isempty(point.candidates) ? :no_candidate : :success,
        candidates=point.candidates,
        diagnostics=point.diagnostics,
        n_candidates=length(point.candidates),
    )
end

function select_operating_point_branch(
    point_result;
    policy::Symbol=:nearest,
    reference_coordinate::Union{Nothing,Real}=nothing,
)
    candidates = point_result.candidates
    isempty(candidates) && return nothing

    idx = if policy == :low
        argmin(getproperty.(candidates, :branch_coordinate))
    elseif policy == :high
        argmax(getproperty.(candidates, :branch_coordinate))
    elseif policy == :nearest
        isnothing(reference_coordinate) && error("reference_coordinate is required for policy=:nearest")
        argmin(abs.(getproperty.(candidates, :branch_coordinate) .- Float64(reference_coordinate)))
    else
        error("unsupported branch-selection policy=$(policy)")
    end

    candidate = candidates[idx]
    return (
        branch_index=idx,
        branch_coordinate=candidate.branch_coordinate,
        mdot=candidate.mdot,
        ht_out=candidate.ht_out,
        tau=candidate.tau,
        PR=candidate.pressure_ratio,
        eta=candidate.efficiency,
        power=candidate.power,
        pressure_ratio=candidate.pressure_ratio,
        efficiency=candidate.efficiency,
        residuals=candidate.residuals,
        machine_payload=candidate.machine_payload,
    )
end
