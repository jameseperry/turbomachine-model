using ....Utility: bracket_bisect_roots
using ...Fluids

@inline _primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x

@inline function _positive_finite(x::Real)
    x_f = _primal_value(x)
    return isfinite(x_f) && x_f > 0.0
end

function _station_radius(model::AxialMachineModel, streamtube_radii::AbstractVector{<:Real}, station_index::Int)
    n_rows = length(model.rows)
    1 <= station_index <= (n_rows + 1) || error("station_index out of bounds")
    if station_index == 1
        return Float64(streamtube_radii[1])
    elseif station_index == n_rows + 1
        return Float64(streamtube_radii[end])
    end
    return 0.5 * (Float64(streamtube_radii[station_index - 1]) + Float64(streamtube_radii[station_index]))
end

function _safe_entropy(
    eos::Fluids.AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    (_positive_finite(pressure) && _positive_finite(enthalpy)) || return NaN
    try
        return Fluids.entropy(eos, pressure, enthalpy)
    catch
        return NaN
    end
end

function _safe_pressure_from_enthalpy_entropy(
    eos::Fluids.AbstractEOS,
    enthalpy::Real,
    entropy::Real,
)
    (_positive_finite(enthalpy) && isfinite(_primal_value(entropy))) || return NaN
    try
        pressure = Fluids.pressure_from_enthalpy_entropy(eos, enthalpy, entropy)
        return _positive_finite(pressure) ? Float64(pressure) : NaN
    catch
        return NaN
    end
end

function _safe_temperature(
    eos::Fluids.AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    (_positive_finite(pressure) && _positive_finite(enthalpy)) || return NaN
    try
        temperature = Fluids.temperature(eos, pressure, enthalpy)
        return _positive_finite(temperature) ? Float64(temperature) : NaN
    catch
        return NaN
    end
end

function _safe_density(
    eos::Fluids.AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    (_positive_finite(pressure) && _positive_finite(enthalpy)) || return NaN
    try
        density = Fluids.density(eos, pressure, enthalpy)
        return _positive_finite(density) ? Float64(density) : NaN
    catch
        return NaN
    end
end

function _safe_speed_of_sound(
    eos::Fluids.AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    (_positive_finite(pressure) && _positive_finite(enthalpy)) || return NaN
    try
        a = Fluids.speed_of_sound(eos, pressure, enthalpy)
        return _positive_finite(a) ? Float64(a) : NaN
    catch
        return NaN
    end
end

function _safe_cp(
    eos::Fluids.AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    (_positive_finite(pressure) && _positive_finite(enthalpy)) || return NaN
    try
        cp = Fluids.heat_capacity_cp(eos, pressure, enthalpy)
        return _positive_finite(cp) ? Float64(cp) : NaN
    catch
        return NaN
    end
end

function _thermo_efficiency(
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    pt_out::Float64,
    ht_out::Float64,
)
    (_positive_finite(pt_in) && _positive_finite(ht_in) && _positive_finite(pt_out) && _positive_finite(ht_out)) || return NaN
    h2s = try
        Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    catch
        NaN
    end
    isfinite(h2s) || return NaN
    if h2s >= ht_in
        denom = ht_out - ht_in
        abs(denom) > 1e-12 || return NaN
        return (h2s - ht_in) / denom
    end
    denom = ht_in - h2s
    abs(denom) > 1e-12 || return NaN
    return (ht_in - ht_out) / denom
end

function _station_static_state(
    eos::Fluids.AbstractEOS,
    pt_t::Float64,
    ht_t::Float64,
    s_t::Float64,
    Vx::Float64,
    Vtheta::Float64,
    area::Float64,
    mdot::Float64,
)
    Vmag2 = Vx^2 + Vtheta^2
    h = Fluids.static_enthalpy_from_total(ht_t, sqrt(Vmag2))
    p = _safe_pressure_from_enthalpy_entropy(eos, h, s_t)
    T = _safe_temperature(eos, p, h)
    rho = _safe_density(eos, p, h)
    a = _safe_speed_of_sound(eos, p, h)
    Mach = (_positive_finite(a) && isfinite(Vmag2)) ? sqrt(Vmag2) / a : NaN
    mdot_station = (_positive_finite(rho) && isfinite(Vx) && isfinite(area)) ? rho * Vx * area : NaN
    return (
        h=h,
        p=p,
        T=T,
        rho=rho,
        a=a,
        Mach=Mach,
        V=sqrt(Vmag2),
        mdot_station=mdot_station,
    )
end

struct StreamtubeStationState
    radius::Float64
    area::Float64
    pt_t::Float64
    ht_t::Float64
    s_t::Float64
    p::Float64
    h::Float64
    Tt::Float64
    T::Float64
    rho::Float64
    a::Float64
    Vx::Float64
    Vtheta::Float64
    V::Float64
    Mach::Float64
    mdot_station::Float64
    valid::Bool
    stall::Bool
    choke::Bool
end

struct StreamtubeSolveResult
    pressure_ratio::Float64
    efficiency::Float64
    thermo_efficiency::Float64
    valid::Bool
    stall::Bool
    choke::Bool
    stations::Vector{StreamtubeStationState}
    row_data::Vector{NamedTuple}
    inputs::NamedTuple
end

@inline _station_field(stations::Vector{StreamtubeStationState}, field::Symbol) = [getfield(st, field) for st in stations]

function Base.getproperty(result::StreamtubeSolveResult, name::Symbol)
    if name in fieldnames(StreamtubeSolveResult)
        return getfield(result, name)
    elseif name === :PR
        return getfield(result, :pressure_ratio)
    elseif name === :eta
        return getfield(result, :efficiency)
    elseif name === :pt_t
        return _station_field(getfield(result, :stations), :pt_t)
    elseif name === :ht_t
        return _station_field(getfield(result, :stations), :ht_t)
    elseif name === :Vx
        return _station_field(getfield(result, :stations), :Vx)
    elseif name === :Vtheta
        return _station_field(getfield(result, :stations), :Vtheta)
    elseif name === :station_data
        return [
            (
                station_index=i,
                radius=st.radius,
                area=st.area,
                pt_t=st.pt_t,
                ht_t=st.ht_t,
                s_t=st.s_t,
                p=st.p,
                h=st.h,
                Tt=st.Tt,
                T=st.T,
                rho=st.rho,
                a=st.a,
                Vx=st.Vx,
                Vtheta=st.Vtheta,
                V=st.V,
                Mach=st.Mach,
                mdot_station=st.mdot_station,
                mdot_error=st.mdot_station - getfield(result, :inputs).mdot,
                valid=st.valid,
                stall=st.stall,
                choke=st.choke,
            ) for (i, st) in pairs(getfield(result, :stations))
        ]
    elseif name === :stall_row
        return [st.stall for st in getfield(result, :stations)[2:end]]
    elseif name === :choke_row
        return [st.choke for st in getfield(result, :stations)[2:end]]
    elseif name === :valid_row
        return [st.valid for st in getfield(result, :stations)[2:end]]
    elseif name === :summary
        return (
            pressure_ratio=getfield(result, :pressure_ratio),
            efficiency=getfield(result, :efficiency),
            thermo_efficiency=getfield(result, :thermo_efficiency),
            valid=getfield(result, :valid),
            stall=getfield(result, :stall),
            choke=getfield(result, :choke),
        )
    end
    return getfield(result, name)
end

function Base.propertynames(::StreamtubeSolveResult, private::Bool=false)
    names = (
        fieldnames(StreamtubeSolveResult)...,
        :PR, :eta, :pt_t, :ht_t, :Vx, :Vtheta, :station_data,
        :stall_row, :choke_row, :valid_row, :summary,
    )
    return private ? names : names
end

@inline function _nan_station_state()
    return StreamtubeStationState(
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
        NaN, NaN, NaN, NaN, NaN, NaN,
        false, false, false,
    )
end

function _build_station_state(
    eos::Fluids.AbstractEOS,
    radius::Float64,
    area::Float64,
    pt_t::Float64,
    ht_t::Float64,
    Vx::Float64,
    Vtheta::Float64,
    mdot::Float64,
    valid::Bool=false,
    stall::Bool=false,
    choke::Bool=false,
)
    s_t = _safe_entropy(eos, pt_t, ht_t)
    static = _station_static_state(eos, pt_t, ht_t, s_t, Vx, Vtheta, area, mdot)
    return StreamtubeStationState(
        radius,
        area,
        pt_t,
        ht_t,
        s_t,
        static.p,
        static.h,
        _safe_temperature(eos, pt_t, ht_t),
        static.T,
        static.rho,
        static.a,
        Vx,
        Vtheta,
        static.V,
        static.Mach,
        static.mdot_station,
        valid,
        stall,
        choke,
    )
end

@inline function _with_station_flags(
    state::StreamtubeStationState;
    valid::Bool=state.valid,
    stall::Bool=state.stall,
    choke::Bool=state.choke,
)
    return StreamtubeStationState(
        state.radius,
        state.area,
        state.pt_t,
        state.ht_t,
        state.s_t,
        state.p,
        state.h,
        state.Tt,
        state.T,
        state.rho,
        state.a,
        state.Vx,
        state.Vtheta,
        state.V,
        state.Mach,
        state.mdot_station,
        valid,
        stall,
        choke,
    )
end

function _station_continuity_residual(
    eos::Fluids.AbstractEOS,
    mdot::Float64,
    pt_t::Float64,
    ht_t::Float64,
    s_t::Float64,
    area::Float64,
    Vtheta::Float64,
    Vx::Float64,
)
    Vx > 0.0 || return NaN
    static = _station_static_state(eos, pt_t, ht_t, s_t, Vx, Vtheta, area, mdot)
    _positive_finite(static.rho) || return NaN
    return mdot - static.rho * Vx * area
end

function _row_outlet_state_from_Vx(
    eos::Fluids.AbstractEOS,
    mdot::Float64,
    radius_out::Float64,
    ht_t_in::Float64,
    Vtheta_in::Float64,
    U::Float64,
    k_theta_exit::Float64,
    s_t_row_in::Float64,
    delta_s_t::Float64,
    area_out::Float64,
    Vx_out::Float64,
)
    Vtheta_out = U + Vx_out * k_theta_exit
    ht_t_out = ht_t_in + U * (Vtheta_out - Vtheta_in)
    s_t_out = s_t_row_in + delta_s_t
    pt_t_out = _safe_pressure_from_enthalpy_entropy(eos, ht_t_out, s_t_out)
    state = _build_station_state(
        eos,
        radius_out,
        area_out,
        pt_t_out,
        ht_t_out,
        Vx_out,
        Vtheta_out,
        mdot,
    )
    return StreamtubeStationState(
        state.radius,
        state.area,
        state.pt_t,
        state.ht_t,
        s_t_out,
        state.p,
        state.h,
        state.Tt,
        state.T,
        state.rho,
        state.a,
        state.Vx,
        state.Vtheta,
        state.V,
        state.Mach,
        state.mdot_station,
        state.valid,
        state.stall,
        state.choke,
    )
end

function _is_physically_admissible_outlet(outlet::StreamtubeStationState)
    return (
        _positive_finite(outlet.Vx) &&
        isfinite(outlet.Vtheta) &&
        _positive_finite(outlet.ht_t) &&
        _positive_finite(outlet.pt_t) &&
        _positive_finite(outlet.h) &&
        _positive_finite(outlet.p) &&
        _positive_finite(outlet.rho) &&
        _positive_finite(outlet.a) &&
        isfinite(outlet.Mach) &&
        isfinite(outlet.V) &&
        _positive_finite(outlet.mdot_station)
    )
end

function _select_outlet_candidate(candidates::Vector{<:StreamtubeStationState}, prefer_root::Symbol)
    isempty(candidates) && return nothing
    return prefer_root == :high ? last(candidates) : first(candidates)
end

function _velocity_upper_bound(
    ht_t::Float64,
    Vtheta::Float64,
)
    kinetic_margin = 2 * ht_t - Vtheta^2
    kinetic_margin > 0.0 || return NaN
    return sqrt(kinetic_margin) * (1 - 1e-10)
end

function _select_root(roots::Vector{Float64}, prefer_root::Symbol)
    isempty(roots) && return NaN
    prefer_root == :high && return last(roots)
    return first(roots)
end

function _solve_station_Vx(
    eos::Fluids.AbstractEOS,
    mdot::Float64,
    pt_t::Float64,
    ht_t::Float64,
    s_t::Float64,
    area::Float64,
    Vtheta::Float64;
    prefer::Symbol=:low,
    prior_roots::AbstractVector{<:Real}=Float64[],
)
    x_hi = _velocity_upper_bound(ht_t, Vtheta)
    isfinite(x_hi) || return (converged=false, Vx=NaN, roots=Float64[])
    f = Vx -> _station_continuity_residual(eos, mdot, pt_t, ht_t, s_t, area, Vtheta, Vx)
    roots = bracket_bisect_roots(
        f,
        (1e-10, x_hi);
        n_scan=401,
        root_tol=1e-9,
        max_bisect_iters=80,
        prior_roots=prior_roots,
        dedupe_atol=1e-8,
    )
    isempty(roots) && return (converged=false, Vx=NaN, roots=roots)
    return (converged=true, Vx=_select_root(roots, prefer), roots=roots)
end

function _solve_row_Vx_roots(
    eos::Fluids.AbstractEOS,
    mdot::Float64,
    radius_out::Float64,
    ht_t_in::Float64,
    Vtheta_in::Float64,
    U::Float64,
    k_theta_exit::Float64,
    s_t_row_in::Float64,
    delta_s_t::Float64,
    area_out::Float64,
)
    # The row solve only makes sense while enough total enthalpy remains to
    # support the requested outlet kinetic energy. This defines the root-search
    # interval for the outlet axial velocity.
    x_hi = _velocity_upper_bound(ht_t_in, U)
    isfinite(x_hi) || return Float64[]

    # Solve continuity for the outlet axial velocity using the full trial outlet
    # state reconstructed from each candidate Vx_out.
    function row_residual(Vx_out)
        outlet = _row_outlet_state_from_Vx(
            eos,
            mdot,
            radius_out,
            ht_t_in,
            Vtheta_in,
            U,
            k_theta_exit,
            s_t_row_in,
            delta_s_t,
            area_out,
            Vx_out,
        )
        return _positive_finite(outlet.mdot_station) ? (mdot - outlet.mdot_station) : NaN
    end

    # Multiple roots correspond to multiple admissible outlet-velocity branches
    # for the same row inlet state.
    return bracket_bisect_roots(
        row_residual,
        (1e-10, x_hi);
        n_scan=401,
        root_tol=1e-9,
        max_bisect_iters=80,
        dedupe_atol=1e-8,
    )
end

function _invalid_streamtube_result(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    streamtube_radii::AbstractVector{<:Real};
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    mdot::Float64,
    Vx_inlet::Float64,
    Vtheta_inlet::Float64,
    stall::Bool,
    choke::Bool,
)
    n_rows = length(model.rows)
    n_stations = n_rows + 1
    stations = [_nan_station_state() for _ in 1:n_stations]
    row_data = NamedTuple[]
    return StreamtubeSolveResult(
        NaN,
        NaN,
        NaN,
        false,
        stall,
        choke,
        stations,
        row_data,
        (
            pt_in=pt_in,
            ht_in=ht_in,
            omega=omega,
            mdot=mdot,
            Vx_inlet=Vx_inlet,
            Vtheta_inlet=Vtheta_inlet,
            streamtube_radii=Float64.(streamtube_radii),
        ),
    )
end

function _build_row_data(
    eos::Fluids.AbstractEOS,
    row::AxialRow,
    row_index::Int,
    row_radius::Float64,
    omega::Float64,
    aero_out::NamedTuple,
    st_in::StreamtubeStationState,
    st_out::StreamtubeStationState,
)
    U = row.speed_ratio_to_ref * omega * row_radius
    Wtheta_in = st_in.Vtheta - U
    Wtheta_out = st_out.Vtheta - U
    Win = sqrt(st_in.Vx^2 + Wtheta_in^2)
    Wout = sqrt(st_out.Vx^2 + Wtheta_out^2)
    delta_ht = st_out.ht_t - st_in.ht_t
    k_theta_exit = aero_out.valid ? tan(Float64(aero_out.theta_out)) : NaN
    delta_s_t = aero_out.valid ? Float64(aero_out.delta_s_t) : NaN
    row_thermo_efficiency = abs(U) > 1e-9 ?
        _thermo_efficiency(eos, st_in.pt_t, st_in.ht_t, st_out.pt_t, st_out.ht_t) : NaN
    q_in = (_positive_finite(st_in.rho) && isfinite(st_in.V)) ? 0.5 * st_in.rho * st_in.V^2 : NaN
    stator_loss_coefficient = (abs(U) <= 1e-9 && isfinite(q_in) && q_in > 1e-9) ?
        (st_in.pt_t - st_out.pt_t) / q_in : NaN
    a_in = st_in.a
    a_out = st_out.a
    return (
        row_index=row_index,
        station_in=row_index,
        station_out=row_index + 1,
        r_hub=row.r_hub,
        r_tip=row.r_tip,
        row_radius=row_radius,
        row_annulus_area=row_annulus_area(row),
        area_in=st_in.area,
        area_out=st_out.area,
        theta_metal_in=row.theta_metal_in,
        theta_metal_out=row.theta_metal_out,
        speed_ratio_to_ref=row.speed_ratio_to_ref,
        omega_row=row.speed_ratio_to_ref * omega,
        U=U,
        Vx_in=st_in.Vx,
        Vx_out=st_out.Vx,
        Vtheta_in=st_in.Vtheta,
        Vtheta_out=st_out.Vtheta,
        delta_Vtheta=st_out.Vtheta - st_in.Vtheta,
        Wtheta_in=Wtheta_in,
        Wtheta_out=Wtheta_out,
        W_in=Win,
        W_out=Wout,
        WMach_in=(_positive_finite(a_in) ? Win / a_in : NaN),
        WMach_out=(_positive_finite(a_out) ? Wout / a_out : NaN),
        alpha_in=atan(st_in.Vtheta, st_in.Vx),
        alpha_out=atan(st_out.Vtheta, st_out.Vx),
        beta_in=atan(Wtheta_in, st_in.Vx),
        beta_out=atan(Wtheta_out, st_out.Vx),
        pt_t_in=st_in.pt_t,
        pt_t_out=st_out.pt_t,
        pt_t_ratio=st_out.pt_t / st_in.pt_t,
        ht_t_in=st_in.ht_t,
        ht_t_out=st_out.ht_t,
        Tt_ratio=st_out.Tt / st_in.Tt,
        delta_ht=delta_ht,
        euler_work=U * (st_out.Vtheta - st_in.Vtheta),
        energy_balance_error=delta_ht - U * (st_out.Vtheta - st_in.Vtheta),
        h_in=st_in.h,
        h_out=st_out.h,
        delta_h=st_out.h - st_in.h,
        p_in=st_in.p,
        p_out=st_out.p,
        p_ratio_row=st_out.p / st_in.p,
        Mach_in=st_in.Mach,
        Mach_out=st_out.Mach,
        thermo_efficiency=row_thermo_efficiency,
        psi_row=abs(U) > 1e-9 ? delta_ht / (U^2) : NaN,
        stator_loss_coefficient=stator_loss_coefficient,
        k_theta_exit=k_theta_exit,
        delta_s_t=delta_s_t,
        stall_margin=aero_out.stall_margin,
        incidence=aero_out.diagnostics.incidence,
        deviation=aero_out.diagnostics.deviation,
        theta_in=aero_out.diagnostics.theta_in,
        theta_out=aero_out.diagnostics.theta_out,
        valid=st_out.valid,
        stall=st_out.stall,
        choke=st_out.choke,
    )
end

function _advance_row_dimensional(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    row::AxialRow,
    streamtube_radii::AbstractVector{<:Real},
    row_radius::Float64,
    station_index::Int,
    inlet_state::StreamtubeStationState,
    omega::Float64,
    mdot::Float64,
)
    # Evaluate the blade-row closure at the current inlet state. This supplies
    # the exit flow angle and the modeled stagnation-entropy rise used to
    # advance the downstream total state.
    U = row.speed_ratio_to_ref * omega * row_radius
    aero_out = blade_aero(
        row.aero,
        row.theta_metal_in,
        row.theta_metal_out,
        inlet_state.Vx,
        inlet_state.Vtheta,
        U,
    )
    stall = (aero_out.stall_margin <= 0) || !aero_out.valid
    if !aero_out.valid
        # The row model itself rejected the inlet condition, so return an invalid
        # row result immediately while preserving the diagnostic payload.
        return (
            aero_out=aero_out,
            choke=false,
            outlet_candidates=StreamtubeStationState[],
        )
    end

    if !isfinite(inlet_state.s_t)
        return (
            aero_out=aero_out,
            choke=false,
            outlet_candidates=StreamtubeStationState[],
        )
    end

    k_theta_exit = tan(Float64(aero_out.theta_out))
    delta_s_t = Float64(aero_out.delta_s_t)
    roots = _solve_row_Vx_roots(
        eos,
        mdot,
        _station_radius(model, streamtube_radii, station_index + 1),
        inlet_state.ht_t,
        inlet_state.Vtheta,
        U,
        k_theta_exit,
        inlet_state.s_t,
        delta_s_t,
        station_area(model, station_index + 1),
    )
    if isempty(roots)
        return (
            aero_out=aero_out,
            choke=true,
            outlet_candidates=StreamtubeStationState[],
        )
    end

    # Reconstruct the full outlet state for every root, then filter out roots
    # that violate hard physical admissibility constraints. Branch selection is
    # handled one level up in `streamtube_solve`.
    area_out = station_area(model, station_index + 1)
    radius_out = _station_radius(model, streamtube_radii, station_index + 1)
    outlet_candidates = StreamtubeStationState[]
    for Vx_out in roots
        outlet = _row_outlet_state_from_Vx(
            eos,
            mdot,
            radius_out,
            inlet_state.ht_t,
            inlet_state.Vtheta,
            U,
            k_theta_exit,
            inlet_state.s_t,
            delta_s_t,
            area_out,
            Vx_out,
        )
        _is_physically_admissible_outlet(outlet) || continue
        push!(outlet_candidates, outlet)
    end

    return (
        aero_out=aero_out,
        choke=false,
        outlet_candidates=outlet_candidates,
    )
end

"""
    streamtube_solve(model, eos, streamtube_radii, pt_in, ht_in, omega, Vx_inlet, Vtheta_inlet; prefer_root=:low)

March the axial machine in physical units.

State per station:
- `pt_t`: total pressure
- `ht_t`: total enthalpy
- `Vx`: axial velocity
- `Vtheta`: tangential velocity

The scalar implicit solve remains, with `Vx` as the unknown at each downstream
station.
"""
function streamtube_solve(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    streamtube_radii::AbstractVector{<:Real},
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    Vx_inlet::Real,
    Vtheta_inlet::Real;
    prefer_root::Symbol=:low,
)
    pt_in_f = Float64(pt_in)
    ht_in_f = Float64(ht_in)
    omega_f = Float64(omega)
    Vx_inlet_f = Float64(Vx_inlet)
    Vtheta_inlet_f = Float64(Vtheta_inlet)
    n_rows = length(model.rows)
    n_stations = n_rows + 1

    # Validate the requested streamtube geometry before doing any thermodynamics.
    # Each row solve assumes the specified streamtube radius lies inside that row's
    # annulus, because blade speed and station area are evaluated at that radius.
    length(streamtube_radii) == n_rows || error("streamtube_radii length must match number of rows")
    radii = Float64.(streamtube_radii)
    for (k, row) in pairs(model.rows)
        row.r_hub <= radii[k] <= row.r_tip ||
            error("streamtube_radii[$k]=$(radii[k]) must lie in [r_hub, r_tip]=[$(row.r_hub), $(row.r_tip)]")
    end

    # Reject obviously invalid inlet conditions up front and return a shaped
    # diagnostic result rather than letting the EOS or row march fail deeper in
    # the solve.
    (_positive_finite(pt_in_f) && _positive_finite(ht_in_f) && _positive_finite(Vx_inlet_f)) ||
        return _invalid_streamtube_result(model, eos, radii; pt_in=pt_in_f, ht_in=ht_in_f, omega=omega_f, mdot=NaN, Vx_inlet=Vx_inlet_f, Vtheta_inlet=Vtheta_inlet_f, stall=true, choke=false)

    # Allocate the physical per-station state. The dimensional march advances a
    # typed station object and the legacy arrays are derived afterward.
    stations = [_nan_station_state() for _ in 1:n_stations]

    # Seed station 1 directly from the specified inlet total state and velocity
    # vector. Mass flow becomes a derived quantity of the inlet station state.
    inlet = _build_station_state(
        eos,
        _station_radius(model, radii, 1),
        station_area(model, 1),
        pt_in_f,
        ht_in_f,
        Vx_inlet_f,
        Vtheta_inlet_f,
        NaN,
        true,
        false,
        false,
    )
    mdot_f = inlet.mdot_station
    _positive_finite(mdot_f) ||
        return _invalid_streamtube_result(model, eos, radii; pt_in=pt_in_f, ht_in=ht_in_f, omega=omega_f, mdot=mdot_f, Vx_inlet=Vx_inlet_f, Vtheta_inlet=Vtheta_inlet_f, stall=false, choke=true)
    stations[1] = inlet

    # Cache row aero diagnostics separately. The row march may terminate early on
    # invalid/choked rows, but the result still needs complete per-row diagnostics
    # for downstream tooling.
    aero_cache = Vector{NamedTuple}(undef, n_rows)

    # March row-by-row from inlet to outlet. Each row:
    # 1. evaluates the blade aero closure at the current inlet state,
    # 2. uses Euler work to relate Vtheta_out and ht_t_out,
    # 3. solves a scalar continuity equation for Vx_out,
    # 4. updates the downstream station total state.
    for k in 1:n_rows
        row = model.rows[k]
        row_step = _advance_row_dimensional(
            model,
            eos,
            row,
            radii,
            radii[k],
            k,
            stations[k],
            omega_f,
            mdot_f,
        )
        aero_cache[k] = row_step.aero_out

        selected_outlet = _select_outlet_candidate(row_step.outlet_candidates, prefer_root)
        if isnothing(selected_outlet)
            break
        end

        stations[k + 1] = _with_station_flags(
            selected_outlet;
            valid=true,
            stall=(row_step.aero_out.stall_margin <= 0) || !row_step.aero_out.valid,
            choke=row_step.choke,
        )
    end

    # Fill any rows skipped after an early termination with placeholder aero
    # diagnostics so consumers can still render a complete machine layout.
    for k in 1:n_rows
        if !isassigned(aero_cache, k)
            aero_cache[k] = (
                k_theta_exit=NaN,
                delta_s_t=NaN,
                stall_margin=NaN,
                valid=false,
                diagnostics=(
                    incidence=NaN,
                    deviation=NaN,
                    theta_in=NaN,
                    theta_out=NaN,
                    theta_metal_in=model.rows[k].theta_metal_in,
                    theta_metal_out=model.rows[k].theta_metal_out,
                ),
            )
        end
    end

    # Reconstruct full station diagnostics directly from the dimensional marched
    # state, then derive row diagnostics from adjacent stations plus the cached
    # aero / row-core solve data.
    pt_t = [st.pt_t for st in stations]
    ht_t = [st.ht_t for st in stations]
    Vtheta = [st.Vtheta for st in stations]
    Vx = [st.Vx for st in stations]
    row_data = Vector{NamedTuple}(undef, n_rows)
    for k in 1:n_rows
        row_data[k] = _build_row_data(
            eos,
            model.rows[k],
            k,
            radii[k],
            omega_f,
            aero_cache[k],
            stations[k],
            stations[k + 1],
        )
    end

    # Collapse the row-wise diagnostics to machine-level validity and performance
    # metrics. Efficiency is now taken directly from the EOS-based thermo
    # calculation rather than a separate nondimensional proxy.
    model_valid = all(st.valid for st in stations[2:end]) && all(isfinite.(pt_t)) && all(isfinite.(ht_t)) && all(isfinite.(Vx)) && all(isfinite.(Vtheta))
    model_stall = any(st.stall for st in stations[2:end])
    model_choke = any(st.choke for st in stations[2:end])
    pressure_ratio = pt_t[end] / pt_t[1]
    thermo_eta = _thermo_efficiency(eos, pt_t[1], ht_t[1], pt_t[end], ht_t[end])
    eta = thermo_eta
    return StreamtubeSolveResult(
        pressure_ratio,
        eta,
        thermo_eta,
        model_valid && isfinite(pressure_ratio) && pressure_ratio > 0 && isfinite(eta),
        model_stall,
        model_choke,
        stations,
        row_data,
        (
            pt_in=pt_in_f,
            ht_in=ht_in_f,
            omega=omega_f,
            mdot=mdot_f,
            Vx_inlet=Vx_inlet_f,
            Vtheta_inlet=Vtheta_inlet_f,
            streamtube_radii=radii,
        ),
    )
end

function streamtube_solve_from_mdot(
    model::AxialMachineModel,
    eos::Fluids.AbstractEOS,
    streamtube_radii::AbstractVector{<:Real},
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    mdot::Real,
    Vtheta_inlet::Real;
    prefer_root::Symbol=:low,
)
    pt_in_f = Float64(pt_in)
    ht_in_f = Float64(ht_in)
    omega_f = Float64(omega)
    mdot_f = Float64(mdot)
    Vtheta_inlet_f = Float64(Vtheta_inlet)
    radii = Float64.(streamtube_radii)

    (_positive_finite(pt_in_f) && _positive_finite(ht_in_f) && _positive_finite(mdot_f)) ||
        return _invalid_streamtube_result(model, eos, radii; pt_in=pt_in_f, ht_in=ht_in_f, omega=omega_f, mdot=mdot_f, Vx_inlet=NaN, Vtheta_inlet=Vtheta_inlet_f, stall=true, choke=false)

    inlet = _solve_station_Vx(
        eos,
        mdot_f,
        pt_in_f,
        ht_in_f,
        _safe_entropy(eos, pt_in_f, ht_in_f),
        station_area(model, 1),
        Vtheta_inlet_f;
        prefer=prefer_root,
    )
    inlet.converged ||
        return _invalid_streamtube_result(model, eos, radii; pt_in=pt_in_f, ht_in=ht_in_f, omega=omega_f, mdot=mdot_f, Vx_inlet=NaN, Vtheta_inlet=Vtheta_inlet_f, stall=false, choke=true)

    return streamtube_solve(
        model,
        eos,
        radii,
        pt_in_f,
        ht_in_f,
        omega_f,
        inlet.Vx,
        Vtheta_inlet_f;
        prefer_root=prefer_root,
    )
end
