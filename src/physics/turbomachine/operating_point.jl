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

function _operating_point_mdot_bounds(
    map::TabulatedPerformanceMap,
    omega::Float64,
    Tt_in::Float64,
    pt_in::Float64,
)
    domain = performance_map_domain(map, Tt_in, pt_in)
    mdot_lo = domain.mdot_flow_range.min(omega)
    mdot_hi = domain.mdot_flow_range.max(omega)
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

function _operating_point_streamtube_inputs(
    model::AxialMachine.AxialMachineModel;
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    mdot::Float64,
    streamtube_radii::AbstractVector{<:Real},
    nu_theta_inlet::Float64,
)
    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    a0_in = sqrt(model.gamma * model.gas_constant * Tt_in)
    rho0_in = pt_in / (model.gas_constant * Tt_in)
    area_in = AxialMachine.station_area(model, 1)
    nu_x_inlet = mdot / (rho0_in * area_in * a0_in)
    m_tip = omega * model.r_tip_ref / a0_in
    return (
        Tt_in=Tt_in,
        a0_in=a0_in,
        rho0_in=rho0_in,
        area_in=area_in,
        m_tip=m_tip,
        nu_x_inlet=nu_x_inlet,
        nu_theta_inlet=nu_theta_inlet,
        streamtube_radii=Float64.(streamtube_radii),
    )
end

function _reconstruct_station_physical_data(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS,
    station_data::Vector{<:NamedTuple};
    pt_in::Float64,
    Tt_in::Float64,
    a0_in::Float64,
    mdot::Float64,
)
    gamma = model.gamma
    return [
        begin
            Vx = st.nu_x * a0_in
            Vtheta = st.nu_theta * a0_in
            Vmag2 = Vx^2 + Vtheta^2
            Tt = st.tau * Tt_in
            ht = Fluids.enthalpy_from_temperature(eos, Tt)
            pt = st.pi * pt_in
            static_ratio = 1 - ((gamma - 1) / (2 * st.tau)) * (st.nu_x^2 + st.nu_theta^2)
            T = static_ratio > 0 ? Tt * static_ratio : NaN
            p = static_ratio > 0 ? pt * static_ratio^(gamma / (gamma - 1)) : NaN
            a = (isfinite(T) && T > 0) ? Fluids.speed_of_sound_from_temperature(eos, T) : NaN
            Mach = isfinite(a) && a > 0 ? sqrt(Vmag2) / a : NaN
            rho = (isfinite(p) && isfinite(T) && T > 0) ? Fluids.density_from_pressure_temperature(eos, p, T) : NaN
            h = isfinite(ht) ? Fluids.static_enthalpy_from_total(ht, sqrt(Vmag2)) : NaN
            mdot_station = (isfinite(rho) && isfinite(Vx)) ? rho * Vx * st.area : NaN
            (
                station_index=st.station_index,
                radius=st.radius,
                area=st.area,
                tau=st.tau,
                pi=st.pi,
                nu_x=st.nu_x,
                nu_theta=st.nu_theta,
                Vx=Vx,
                Vtheta=Vtheta,
                V=sqrt(Vmag2),
                Tt=Tt,
                pt=pt,
                ht=ht,
                h=h,
                T=T,
                p=p,
                rho=rho,
                Mach=Mach,
                mdot_station=mdot_station,
                mdot_error=mdot_station - mdot,
            )
        end
        for st in station_data
    ]
end

function _reconstruct_row_physical_data(
    eos::Fluids.AbstractEOS,
    row_data::Vector{<:NamedTuple},
    physical_station_data::Vector{<:NamedTuple};
    a0_in::Float64,
    omega::Float64,
)
    return [
        let
            st_in = physical_station_data[row.station_in]
            st_out = physical_station_data[row.station_out]
            omega_row = row.speed_ratio_to_ref * omega
            U = row.nu_u * a0_in
            Wtheta_in = st_in.Vtheta - row.nu_u * a0_in
            Wtheta_out = st_out.Vtheta - row.nu_u * a0_in
            Win = sqrt(st_in.Vx^2 + Wtheta_in^2)
            Wout = sqrt(st_out.Vx^2 + Wtheta_out^2)
            ht_in = st_in.ht
            ht_out = st_out.ht
            delta_ht = ht_out - ht_in
            euler_work = row.nu_u == 0.0 ? 0.0 : U * (st_out.Vtheta - st_in.Vtheta)
            row_thermo_efficiency = abs(euler_work) > 1e-9 ?
                _thermo_efficiency(
                    eos,
                    st_in.pt,
                    ht_in,
                    st_out.pt,
                    ht_out,
                ) : NaN
            h_in = st_in.h
            h_out = st_out.h
            delta_h = h_out - h_in
            pt_ratio_row = st_out.pt / st_in.pt
            Tt_ratio_row = st_out.Tt / st_in.Tt
            p_ratio_row = st_out.p / st_in.p
            a_in = (isfinite(st_in.T) && st_in.T > 0) ? Fluids.speed_of_sound_from_temperature(eos, st_in.T) : NaN
            a_out = (isfinite(st_out.T) && st_out.T > 0) ? Fluids.speed_of_sound_from_temperature(eos, st_out.T) : NaN
            WMach_in = isfinite(a_in) && a_in > 0 ? Win / a_in : NaN
            WMach_out = isfinite(a_out) && a_out > 0 ? Wout / a_out : NaN
            psi_row = abs(U) > 1e-9 ? delta_ht / (U^2) : NaN
            q_in = (isfinite(st_in.rho) && isfinite(st_in.V)) ? 0.5 * st_in.rho * st_in.V^2 : NaN
            stator_loss_coefficient = (abs(U) <= 1e-9 && isfinite(q_in) && q_in > 1e-9) ?
                (st_in.pt - st_out.pt) / q_in : NaN
            (
                row_index=row.row_index,
                station_in=row.station_in,
                station_out=row.station_out,
                r_hub=row.r_hub,
                r_tip=row.r_tip,
                row_radius=row.row_radius,
                row_annulus_area=row.row_annulus_area,
                area_in=row.area_in,
                area_out=row.area_out,
                theta_metal_in=row.theta_metal_in,
                theta_metal_out=row.theta_metal_out,
                speed_ratio_to_ref=row.speed_ratio_to_ref,
                omega_row=omega_row,
                nu_u=row.nu_u,
                U=U,
                nu_x_in=row.nu_x_in,
                nu_x_out=row.nu_x_out,
                nu_theta_in=row.nu_theta_in,
                nu_theta_out=row.nu_theta_out,
                delta_nu_theta=row.delta_nu_theta,
                tau_in=row.tau_in,
                tau_out=row.tau_out,
                delta_tau=row.delta_tau,
                pi_in=row.pi_in,
                pi_out=row.pi_out,
                delta_pi=row.delta_pi,
                k_theta_exit=row.k_theta_exit,
                delta_s_hat=row.delta_s_hat,
                stall_margin=row.stall_margin,
                incidence=row.incidence,
                deviation=row.deviation,
                theta_in=row.theta_in,
                theta_out=row.theta_out,
                Vx_in=st_in.Vx,
                Vx_out=st_out.Vx,
                Vtheta_in=st_in.Vtheta,
                Vtheta_out=st_out.Vtheta,
                Wtheta_in=Wtheta_in,
                Wtheta_out=Wtheta_out,
                W_in=Win,
                W_out=Wout,
                alpha_in=atan(st_in.Vtheta, st_in.Vx),
                alpha_out=atan(st_out.Vtheta, st_out.Vx),
                beta_in=atan(Wtheta_in, st_in.Vx),
                beta_out=atan(Wtheta_out, st_out.Vx),
                ht_in=ht_in,
                ht_out=ht_out,
                delta_ht=delta_ht,
                h_in=h_in,
                h_out=h_out,
                delta_h=delta_h,
                euler_work=euler_work,
                energy_balance_error=delta_ht - euler_work,
                thermo_efficiency=row_thermo_efficiency,
                pt_ratio_row=pt_ratio_row,
                Tt_ratio_row=Tt_ratio_row,
                p_ratio_row=p_ratio_row,
                WMach_in=WMach_in,
                WMach_out=WMach_out,
                psi_row=psi_row,
                stator_loss_coefficient=stator_loss_coefficient,
                pt_in=st_in.pt,
                pt_out=st_out.pt,
                Tt_in=st_in.Tt,
                Tt_out=st_out.Tt,
                p_in=st_in.p,
                p_out=st_out.p,
                T_in=st_in.T,
                T_out=st_out.T,
                Mach_in=st_in.Mach,
                Mach_out=st_out.Mach,
                valid=row.valid,
                stall=row.stall,
                choke=row.choke,
            )
        end
        for row in row_data
    ]
end

function _reconstruct_streamtube_diagnostics(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS,
    raw::NamedTuple,
    inputs::NamedTuple;
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    mdot::Float64,
    operating_point_summary::NamedTuple,
)
    physical_station_data = _reconstruct_station_physical_data(
        model,
        eos,
        raw.station_data;
        pt_in=pt_in,
        Tt_in=inputs.Tt_in,
        a0_in=inputs.a0_in,
        mdot=mdot,
    )
    physical_row_data = _reconstruct_row_physical_data(
        eos,
        raw.row_data,
        physical_station_data;
        a0_in=inputs.a0_in,
        omega=omega,
    )
    return (
        inputs=inputs,
        nondimensional=(
            station_data=raw.station_data,
            row_data=raw.row_data,
            tau=raw.tau,
            pi=raw.pi,
            nu_theta=raw.nu_theta,
            nu_x=raw.nu_x,
        ),
        physical=(
            station_data=physical_station_data,
            row_data=physical_row_data,
        ),
        summary=operating_point_summary,
    )
end

function _stage_diagnostics_from_candidate(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS,
    candidate::NamedTuple;
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    streamtube_radii::AbstractVector{<:Real},
    nu_theta_inlet::Float64,
    prefer_root::Symbol,
)
    inputs = _operating_point_streamtube_inputs(
        model;
        eos=eos,
        pt_in=pt_in,
        ht_in=ht_in,
        omega=omega,
        mdot=Float64(candidate.mdot),
        streamtube_radii=streamtube_radii,
        nu_theta_inlet=nu_theta_inlet,
    )
    raw = AxialMachine.streamtube_solve(
        model,
        inputs.streamtube_radii,
        inputs.m_tip,
        inputs.nu_x_inlet,
        inputs.nu_theta_inlet;
        prefer_root=prefer_root,
    )
    ht_out_streamtube = Fluids.enthalpy_from_temperature(eos, raw.tau[end] * inputs.Tt_in)
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
        raw.PR * pt_in,
        ht_out_streamtube,
    )
    return _reconstruct_streamtube_diagnostics(
        model,
        eos,
        raw,
        inputs;
        pt_in=pt_in,
        ht_in=ht_in,
        omega=omega,
        mdot=Float64(candidate.mdot),
        operating_point_summary=(
            operating_point_pressure_ratio=Float64(candidate.pressure_ratio),
            operating_point_efficiency=Float64(candidate.efficiency),
            operating_point_thermo_efficiency=operating_point_thermo_efficiency,
            operating_point_ht_out=Float64(candidate.ht_out),
            operating_point_tau=Float64(candidate.tau),
            operating_point_power=Float64(candidate.power),
            streamtube_pressure_ratio=raw.PR,
            streamtube_efficiency=raw.eta,
            streamtube_thermo_efficiency=streamtube_thermo_efficiency,
            streamtube_ht_out=ht_out_streamtube,
            streamtube_tau=tau_streamtube,
            streamtube_power=power_streamtube,
            pressure_ratio_error=raw.PR - Float64(candidate.pressure_ratio),
            efficiency_error=raw.eta - Float64(candidate.efficiency),
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
    diagnose_axial_operating_point(model, eos; pt_in, ht_in, omega, mdot, streamtube_radii=meanline_radii(model), nu_theta_inlet=0.0, prefer_root=:low)

Run the axial-machine streamtube solver directly at one physical operating point
and return detailed non-dimensional and reconstructed physical row/station data.
"""
function diagnose_axial_operating_point(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    mdot::Real,
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    nu_theta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
)
    pt_in_f = Float64(pt_in)
    ht_in_f = Float64(ht_in)
    omega_f = Float64(omega)
    mdot_f = Float64(mdot)
    inputs = _operating_point_streamtube_inputs(
        model;
        eos=eos,
        pt_in=pt_in_f,
        ht_in=ht_in_f,
        omega=omega_f,
        mdot=mdot_f,
        streamtube_radii=streamtube_radii,
        nu_theta_inlet=Float64(nu_theta_inlet),
    )
    raw = AxialMachine.streamtube_solve(
        model,
        inputs.streamtube_radii,
        inputs.m_tip,
        inputs.nu_x_inlet,
        inputs.nu_theta_inlet;
        prefer_root=prefer_root,
    )
    Tt_out = raw.tau[end] * inputs.Tt_in
    ht_out = Fluids.enthalpy_from_temperature(eos, Tt_out)
    tau_shaft = mdot_f * (ht_out - ht_in_f) / _safe_omega_divisor(omega_f)
    power = tau_shaft * omega_f
    thermo_efficiency = _thermo_efficiency(eos, pt_in_f, ht_in_f, raw.PR * pt_in_f, ht_out)
    return _reconstruct_streamtube_diagnostics(
        model,
        eos,
        raw,
        inputs;
        pt_in=pt_in_f,
        ht_in=ht_in_f,
        omega=omega_f,
        mdot=mdot_f,
        operating_point_summary=(
            mdot=mdot_f,
            omega=omega_f,
            pressure_ratio=raw.PR,
            efficiency=raw.eta,
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
    replay_operating_point_with_streamtube(model, eos, candidate; pt_in, ht_in, omega, streamtube_radii=meanline_radii(model), nu_theta_inlet=0.0, prefer_root=:low)

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
    nu_theta_inlet::Real=0.0,
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
        nu_theta_inlet=Float64(nu_theta_inlet),
        prefer_root=prefer_root,
    )
end

function replay_operating_point_with_streamtube(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    mdot::Real,
    pressure_ratio::Real,
    efficiency::Real,
    ht_out::Real,
    tau::Real,
    power::Real=(Float64(tau) * Float64(omega)),
    streamtube_radii::AbstractVector{<:Real}=AxialMachine.meanline_radii(model),
    nu_theta_inlet::Real=0.0,
    prefer_root::Symbol=:low,
)
    candidate = (
        branch_coordinate=Float64(mdot),
        mdot=Float64(mdot),
        ht_out=Float64(ht_out),
        tau=Float64(tau),
        pressure_ratio=Float64(pressure_ratio),
        efficiency=Float64(efficiency),
        power=Float64(power),
    )
    return replay_operating_point_with_streamtube(
        model,
        eos,
        candidate;
        pt_in=pt_in,
        ht_in=ht_in,
        omega=omega,
        streamtube_radii=streamtube_radii,
        nu_theta_inlet=nu_theta_inlet,
        prefer_root=prefer_root,
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
