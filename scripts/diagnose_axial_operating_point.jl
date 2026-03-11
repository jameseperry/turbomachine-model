#!/usr/bin/env julia

using ArgParse
using Printf
using TurboMachineModel

const TP = TurboMachineModel.Physics
const TM = TP.Turbomachine
const Axial = TM.AxialMachine
const Fluids = TP.Fluids
const RAD2DEG = 180.0 / pi

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _build_parser()
    settings = ArgParseSettings(
        prog="diagnose_axial_operating_point.jl",
        description="Run one axial-machine operating point through the streamtube solver and print detailed stage diagnostics.",
    )

    @add_arg_table! settings begin
        "model_path"
            help = "input axial model TOML path"
            required = true
        "--model-group"
            help = "input TOML group/table name"
            arg_type = String
            default = "axial_model"
        "--fluid"
            help = "working fluid EOS key"
            arg_type = String
            default = "air"
        "--pt-in"
            help = "inlet total pressure [Pa]"
            arg_type = Float64
            required = true
        "--tt-in"
            help = "inlet total temperature [K]"
            arg_type = Float64
            required = true
        "--omega"
            help = "shaft speed [rad/s]"
            arg_type = Float64
            required = true
        "--vx-inlet"
            help = "inlet axial velocity [m/s]"
            arg_type = Float64
            required = true
        "--vtheta-inlet"
            help = "inlet tangential velocity [m/s]"
            arg_type = Float64
            default = 0.0
        "--prefer-root"
            help = "streamtube root preference when multiple local roots exist: low or high"
            arg_type = String
            default = "low"
        "--output"
            help = "optional output text file path; defaults to stdout"
            arg_type = String
    end

    return settings
end

function _resolve_eos(fluid_key::AbstractString)
    key = Symbol(lowercase(strip(fluid_key)))
    registry = Fluids.ideal_EOS()
    haskey(registry, key) || error("unsupported fluid=$(fluid_key)")
    return registry[key]
end

function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    x isa Integer && return string(x)
    x isa AbstractFloat || return string(x)
    if isnan(x)
        return "NaN"
    elseif isinf(x)
        return x > 0 ? "Inf" : "-Inf"
    elseif abs(x) >= 1e4 || (abs(x) > 0 && abs(x) < 1e-3)
        return @sprintf("%.4e", x)
    else
        return @sprintf("%.6f", x)
    end
end

function _print_named_tuple_block(io::IO, title::AbstractString, data::NamedTuple)
    println(io, title)
    keys_list = collect(keys(data))
    key_width = maximum(length.(string.(keys_list)))
    for key in keys_list
        println(io, "  ", rpad(String(key), key_width), "  ", _fmt(getproperty(data, key)))
    end
    println(io)
end

function _print_table(io::IO, title::AbstractString, rows::AbstractVector{<:NamedTuple}, columns::Vector{Symbol})
    println(io, title)
    isempty(rows) && (println(io, "  <empty>\n"); return)
    widths = Dict{Symbol,Int}()
    for col in columns
        widths[col] = maximum(vcat(length(String(col)), [length(_fmt(getproperty(row, col))) for row in rows]))
    end
    println(io, "  ", join([rpad(String(col), widths[col]) for col in columns], "  "))
    for row in rows
        println(io, "  ", join([lpad(_fmt(getproperty(row, col)), widths[col]) for col in columns], "  "))
    end
    println(io)
end

_row_as_named_tuple(row) = (; (name => getproperty(row, name) for name in propertynames(row))...)

function _convert_angle_columns(rows::AbstractVector, angle_cols::Vector{Symbol})
    return [
        begin
            base = _row_as_named_tuple(row)
            updates = (; (Symbol(String(col) * "_deg") => getproperty(row, col) * RAD2DEG for col in angle_cols)...)
            merge(base, updates)
        end
        for row in rows
    ]
end

function _thermo_efficiency_row(eos, pt_in, ht_in, pt_out, ht_out, U)
    abs(U) > 1e-9 || return NaN
    if !(isfinite(pt_in) && isfinite(ht_in) && isfinite(pt_out) && isfinite(ht_out))
        return NaN
    end
    if pt_in <= 0 || pt_out <= 0
        return NaN
    end
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    if pt_out >= pt_in
        denom = ht_out - ht_in
        return abs(denom) > 1e-12 ? (h2s - ht_in) / denom : NaN
    end
    denom = ht_in - h2s
    return abs(denom) > 1e-12 ? (ht_in - ht_out) / denom : NaN
end

function _row_display_rows(model, eos, inputs::NamedTuple, stations, rows)
    radii = inputs.streamtube_radii
    omega = inputs.omega
    return [
        begin
            row_model = model.rows[row.row_index]
            st_in = stations[row.row_index]
            st_out = stations[row.row_index + 1]
            row_radius = Axial.row_streamtube_radius(radii, row.row_index)
            U = Axial.row_angular_speed(row_model, omega) * row_radius
            Wtheta_in = st_in.Vtheta - U
            Wtheta_out = st_out.Vtheta - U
            W_in = sqrt(st_in.Vx^2 + Wtheta_in^2)
            W_out = sqrt(st_out.Vx^2 + Wtheta_out^2)
            q_in = (isfinite(st_in.rho) && st_in.rho > 0 && isfinite(st_in.V)) ? 0.5 * st_in.rho * st_in.V^2 : NaN
            merge(
                _row_as_named_tuple(row),
                (
                    r_hub=row_model.r_hub,
                    r_tip=row_model.r_tip,
                    row_radius=row_radius,
                    omega_row=Axial.row_angular_speed(row_model, omega),
                    U=U,
                    theta_metal_in=row_model.theta_metal_in,
                    theta_metal_out=row_model.theta_metal_out,
                    theta_in=row.aero.theta_in,
                    theta_out=row.aero.theta_out,
                    regime=row.aero.regime,
                    incidence=row.aero.incidence,
                    deviation=row.aero.deviation,
                    delta_s_t=row.aero.delta_s_t,
                    stall_margin=row.aero.stall_margin,
                    n_outlet_candidates=length(row.outlet_candidates),
                    selected_candidate_index=row.selected_candidate_index,
                    valid=st_out.valid,
                    stall=st_out.stall,
                    Vx_in=st_in.Vx,
                    Vx_out=st_out.Vx,
                    Vtheta_in=st_in.Vtheta,
                    Vtheta_out=st_out.Vtheta,
                    W_in=W_in,
                    W_out=W_out,
                    WMach_in=(isfinite(st_in.a) && st_in.a > 0) ? W_in / st_in.a : NaN,
                    WMach_out=(isfinite(st_out.a) && st_out.a > 0) ? W_out / st_out.a : NaN,
                    alpha_in=atan(st_in.Vtheta, st_in.Vx),
                    alpha_out=atan(st_out.Vtheta, st_out.Vx),
                    beta_in=atan(Wtheta_in, st_in.Vx),
                    beta_out=atan(Wtheta_out, st_out.Vx),
                    delta_ht=st_out.ht_t - st_in.ht_t,
                    euler_work=U * (st_out.Vtheta - st_in.Vtheta),
                    thermo_efficiency=_thermo_efficiency_row(eos, st_in.pt_t, st_in.ht_t, st_out.pt_t, st_out.ht_t, U),
                    psi_row=abs(U) > 1e-9 ? (st_out.ht_t - st_in.ht_t) / (U^2) : NaN,
                    pt_t_ratio=st_out.pt_t / st_in.pt_t,
                    Tt_ratio=st_out.Tt / st_in.Tt,
                    p_ratio_row=st_out.p / st_in.p,
                    stator_loss_coefficient=(abs(U) <= 1e-9 && isfinite(q_in) && q_in > 1e-9) ? (st_in.pt_t - st_out.pt_t) / q_in : NaN,
                ),
            )
        end
        for row in rows
    ]
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    model_group = String(something(_parsed_opt(parsed, "model_group", "model-group"), "axial_model"))
    prefer_root = Symbol(lowercase(String(parsed["prefer-root"])))
    prefer_root in (:low, :high) || error("prefer-root must be low or high")

    model = Axial.read_toml(TM.AxialModel, parsed["model_path"]; group=model_group)
    eos = _resolve_eos(parsed["fluid"])
    pt_in = Float64(parsed["pt-in"])
    Tt_in = Float64(parsed["tt-in"])
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)
    omega = Float64(parsed["omega"])
    Vx_inlet = Float64(parsed["vx-inlet"])

    diagnostics = TM.diagnose_axial_operating_point(
        model,
        eos;
        pt_in=pt_in,
        ht_in=ht_in,
        omega=omega,
        Vx_inlet=Vx_inlet,
        Vtheta_inlet=Float64(parsed["vtheta-inlet"]),
        prefer_root=prefer_root,
    )

    output_path = _parsed_opt(parsed, "output", "output")
    io = if isnothing(output_path)
        stdout
    else
        mkpath(dirname(String(output_path)))
        open(String(output_path), "w")
    end

    try
        _print_named_tuple_block(io, "Operating Point Summary", diagnostics.summary)
        _print_named_tuple_block(io, "Derived Streamtube Inputs", diagnostics.inputs)

        _print_table(
        io,
        "Station Diagnostics (Physical)",
        diagnostics.physical.station_data,
        [:station_index, :radius, :area, :Tt, :pt_t, :T, :p, :rho, :ht_t, :h, :Vx, :Vtheta, :V, :Mach, :mdot_station, :mdot_error],
        )
        phys_row_display = _convert_angle_columns(
            _row_display_rows(model, eos, diagnostics.inputs, diagnostics.physical.station_data, diagnostics.physical.row_data),
            [:incidence, :deviation, :theta_metal_in, :theta_metal_out, :theta_in, :theta_out, :alpha_in, :alpha_out, :beta_in, :beta_out],
        )
        _print_table(
        io,
        "Row Diagnostics (Physical)",
        phys_row_display,
        [:row_index, :n_outlet_candidates, :selected_candidate_index, :regime, :r_hub, :r_tip, :row_radius, :omega_row, :U, :theta_metal_in_deg, :theta_metal_out_deg, :theta_in_deg, :theta_out_deg, :incidence_deg, :deviation_deg, :delta_s_t, :stall_margin, :Vx_in, :Vx_out, :Vtheta_in, :Vtheta_out, :W_in, :W_out, :WMach_in, :WMach_out, :alpha_in_deg, :alpha_out_deg, :beta_in_deg, :beta_out_deg, :delta_ht, :euler_work, :psi_row, :thermo_efficiency, :pt_t_ratio, :Tt_ratio, :p_ratio_row, :stator_loss_coefficient, :valid, :stall, :choke],
        )
    finally
        io === stdout || close(io)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
