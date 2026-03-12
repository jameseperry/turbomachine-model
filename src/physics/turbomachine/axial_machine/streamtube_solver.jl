using ....Utility: bracket_bisect_roots
using ...Fluids

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

Base.@kwdef struct StreamtubeRowState
    row_index::Int
    aero::BladeAeroResult
    outlet_candidates::Vector{StreamtubeStationState}
    selected_candidate_index::Int
    choke::Bool
end

struct StreamtubeSolveResult
    pressure_ratio::Float64
    efficiency::Float64
    thermo_efficiency::Float64
    valid::Bool
    stall::Bool
    choke::Bool
    status::Symbol
    stations::Vector{StreamtubeStationState}
    row_data::Vector{StreamtubeRowState}
end

const _W_MACH_CHOKE_LIMIT = 0.98

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
                mdot_error=st.mdot_station - getfield(result, :stations)[1].mdot_station,
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
            status=getfield(result, :status),
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
    return prefer_root == :high ? length(candidates) : 1
end

function _next_row_inlet_WMach(
    row::AxialRow,
    row_radius::Float64,
    inlet_state::StreamtubeStationState,
    omega::Float64,
)
    U = row_angular_speed(row, omega) * row_radius
    Wtheta_in = inlet_state.Vtheta - U
    W_in = hypot(inlet_state.Vx, Wtheta_in)
    return (_positive_finite(inlet_state.a) && isfinite(W_in)) ? (W_in / inlet_state.a) : NaN
end

function _row_outlet_WMach(
    row::AxialRow,
    row_radius::Float64,
    outlet_state::StreamtubeStationState,
    omega::Float64,
)
    U = row_angular_speed(row, omega) * row_radius
    Wtheta_out = outlet_state.Vtheta - U
    W_out = hypot(outlet_state.Vx, Wtheta_out)
    return (_positive_finite(outlet_state.a) && isfinite(W_out)) ? (W_out / outlet_state.a) : NaN
end

function _filter_outlet_candidates_by_next_row_choke(
    model::AxialMachineModel,
    streamtube_radii::AbstractVector{<:Real},
    row_index::Int,
    outlet_candidates::Vector{StreamtubeStationState},
    omega::Float64,
)
    row_index < length(model.rows) || return (outlet_candidates=outlet_candidates, choke=false, status=:normal)
    next_row = model.rows[row_index + 1]
    next_row_radius = row_streamtube_radius(streamtube_radii, row_index + 1)
    filtered = StreamtubeStationState[]
    for candidate in outlet_candidates
        WMa_in = _next_row_inlet_WMach(next_row, next_row_radius, candidate, omega)
        (isfinite(WMa_in) && WMa_in < _W_MACH_CHOKE_LIMIT) || continue
        push!(filtered, candidate)
    end
    return (
        outlet_candidates=filtered,
        choke=isempty(filtered) && !isempty(outlet_candidates),
        status=(isempty(filtered) && !isempty(outlet_candidates)) ? :next_row_inlet_choke : :normal,
    )
end

function _velocity_upper_bound(
    ht_t::Float64,
    Vtheta::Float64,
)
    kinetic_margin = 2 * ht_t - Vtheta^2
    kinetic_margin > 0.0 || return NaN
    return sqrt(kinetic_margin) * (1 - 1e-10)
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
    status::Symbol,
)
    n_rows = length(model.rows)
    n_stations = n_rows + 1
    stations = [_nan_station_state() for _ in 1:n_stations]
    row_data = StreamtubeRowState[]
    return StreamtubeSolveResult(
        NaN,
        NaN,
        NaN,
        false,
        stall,
        choke,
        status,
        stations,
        row_data,
    )
end

_build_row_data(
    row_index::Int,
    aero_out::BladeAeroResult,
    outlet_candidates::Vector{StreamtubeStationState},
    selected_candidate_index::Int,
    choke::Bool,
) =
    StreamtubeRowState(
        row_index=row_index,
        aero=aero_out,
        outlet_candidates=outlet_candidates,
        selected_candidate_index=selected_candidate_index,
        choke=choke,
    )

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
    U = row_angular_speed(row, omega) * row_radius
    aero_out = blade_aero(
        row.aero,
        row.theta_metal_in,
        row.theta_metal_out,
        inlet_state.Vx,
        inlet_state.Vtheta,
        U,
    )
    inlet_WMach = _next_row_inlet_WMach(row, row_radius, inlet_state, omega)
    if !(isfinite(inlet_WMach) && inlet_WMach < _W_MACH_CHOKE_LIMIT)
        return (
            aero_out=aero_out,
            choke=true,
            status=:row_inlet_choke,
            outlet_candidates=StreamtubeStationState[],
        )
    end
    if !isfinite(inlet_state.s_t)
        return (
            aero_out=aero_out,
            choke=false,
            status=:row_invalid_inlet_state,
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
            status=:row_no_roots,
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
        WMach_out = _row_outlet_WMach(row, row_radius, outlet, omega)
        (isfinite(WMach_out) && WMach_out < _W_MACH_CHOKE_LIMIT) || continue
        push!(outlet_candidates, outlet)
    end

    return (
        aero_out=aero_out,
        choke=isempty(outlet_candidates),
        status=isempty(outlet_candidates) ? :row_no_outlet_candidates : :normal,
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
        return _invalid_streamtube_result(model, eos, radii; pt_in=pt_in_f, ht_in=ht_in_f, omega=omega_f, mdot=NaN, Vx_inlet=Vx_inlet_f, Vtheta_inlet=Vtheta_inlet_f, stall=true, choke=false, status=:invalid_inlet)

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
        return _invalid_streamtube_result(model, eos, radii; pt_in=pt_in_f, ht_in=ht_in_f, omega=omega_f, mdot=mdot_f, Vx_inlet=Vx_inlet_f, Vtheta_inlet=Vtheta_inlet_f, stall=false, choke=true, status=:invalid_inlet_massflow)
    stations[1] = inlet

    # Cache row aero diagnostics separately. The row march may terminate early on
    # invalid/choked rows, but the result still needs complete per-row diagnostics
    # for downstream tooling.
    aero_cache = Vector{BladeAeroResult}(undef, n_rows)
    selected_candidate_index_cache = fill(0, n_rows)
    outlet_candidate_cache = [StreamtubeStationState[] for _ in 1:n_rows]
    row_choke_cache = fill(false, n_rows)
    failure_status = :normal

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
        row_candidates = row_step.outlet_candidates
        row_choke = row_step.choke
        if !isempty(row_candidates)
            lookahead = _filter_outlet_candidates_by_next_row_choke(
                model,
                radii,
                k,
                row_candidates,
                omega_f,
            )
            row_candidates = lookahead.outlet_candidates
            row_choke = row_choke || lookahead.choke
            if lookahead.status != :normal
                row_step = (
                    aero_out=row_step.aero_out,
                    choke=row_choke,
                    status=lookahead.status,
                    outlet_candidates=row_candidates,
                )
            end
        end
        outlet_candidate_cache[k] = row_candidates
        row_choke_cache[k] = row_choke

        selected_index = _select_outlet_candidate(row_candidates, prefer_root)
        if isnothing(selected_index)
            failure_status = row_step.status
            break
        end
        selected_candidate_index_cache[k] = selected_index

        stations[k + 1] = _with_station_flags(
            row_candidates[selected_index];
            valid=true,
            stall=(row_step.aero_out.regime == :stall),
            choke=row_choke,
        )
    end

    # Fill any rows skipped after an early termination with placeholder aero
    # diagnostics so consumers can still render a complete machine layout.
    for k in 1:n_rows
        if !isassigned(aero_cache, k)
            aero_cache[k] = BladeAeroResult(
                NaN,
                NaN,
                NaN,
                NaN,
                :stall,
                NaN,
                NaN,
                NaN,
                model.rows[k].theta_metal_in,
                model.rows[k].theta_metal_out,
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
    row_data = Vector{StreamtubeRowState}(undef, n_rows)
    for k in 1:n_rows
        row_data[k] = _build_row_data(
            k,
            aero_cache[k],
            outlet_candidate_cache[k],
            selected_candidate_index_cache[k],
            row_choke_cache[k],
        )
    end

    # Collapse the row-wise diagnostics to machine-level validity and performance
    # metrics. Efficiency is now taken directly from the EOS-based thermo
    # calculation rather than a separate nondimensional proxy.
    model_valid = all(st.valid for st in stations[2:end]) && all(isfinite.(pt_t)) && all(isfinite.(ht_t)) && all(isfinite.(Vx)) && all(isfinite.(Vtheta))
    model_stall = any(st.stall for st in stations[2:end])
    model_choke = any(row_choke_cache)
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
        model_valid ? :normal : failure_status,
        stations,
        row_data,
    )
end
