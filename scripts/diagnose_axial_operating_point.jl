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
        "--mdot"
            help = "mass flow [kg/s]"
            arg_type = Float64
            required = true
        "--nu-theta-inlet"
            help = "inlet non-dimensional swirl velocity"
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

function _convert_angle_columns(rows::AbstractVector{<:NamedTuple}, angle_cols::Vector{Symbol})
    return [
        begin
            updates = (; (Symbol(String(col) * "_deg") => getproperty(row, col) * RAD2DEG for col in angle_cols)...)
            merge(row, updates)
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
    mdot = Float64(parsed["mdot"])

    diagnostics = TM.diagnose_axial_operating_point(
        model,
        eos;
        pt_in=pt_in,
        ht_in=ht_in,
        omega=omega,
        mdot=mdot,
        nu_theta_inlet=Float64(parsed["nu-theta-inlet"]),
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
        "Station Diagnostics (Non-dimensional)",
        diagnostics.nondimensional.station_data,
        [:station_index, :radius, :area, :tau, :pi, :nu_x, :nu_theta],
        )
        _print_table(
        io,
        "Station Diagnostics (Physical)",
        diagnostics.physical.station_data,
        [:station_index, :radius, :area, :Tt, :pt, :T, :p, :rho, :ht, :h, :Vx, :Vtheta, :V, :Mach, :mdot_station, :mdot_error],
        )
        nd_row_display = _convert_angle_columns(
            diagnostics.nondimensional.row_data,
            [:theta_metal_in, :theta_metal_out, :incidence, :deviation, :theta_in, :theta_out],
        )
        _print_table(
        io,
        "Row Diagnostics (Non-dimensional)",
        nd_row_display,
        [:row_index, :row_radius, :row_annulus_area, :theta_metal_in_deg, :theta_metal_out_deg, :nu_u, :nu_x_in, :nu_x_out, :nu_theta_in, :nu_theta_out, :tau_in, :tau_out, :pi_in, :pi_out, :incidence_deg, :deviation_deg, :theta_in_deg, :theta_out_deg, :stall_margin, :valid, :stall, :choke],
        )
        phys_row_display = _convert_angle_columns(
            diagnostics.physical.row_data,
            [:incidence, :deviation, :theta_metal_in, :theta_metal_out, :alpha_in, :alpha_out, :beta_in, :beta_out],
        )
        _print_table(
        io,
        "Row Diagnostics (Physical)",
        phys_row_display,
        [:row_index, :r_hub, :r_tip, :row_radius, :omega_row, :U, :incidence_deg, :deviation_deg, :theta_metal_in_deg, :theta_metal_out_deg, :alpha_in_deg, :alpha_out_deg, :beta_in_deg, :beta_out_deg, :Vx_in, :Vx_out, :Vtheta_in, :Vtheta_out, :W_in, :W_out, :WMach_in, :WMach_out, :ht_in, :ht_out, :delta_ht, :h_in, :h_out, :delta_h, :euler_work, :energy_balance_error, :psi_row, :thermo_efficiency, :pt_ratio_row, :Tt_ratio_row, :p_ratio_row, :stator_loss_coefficient, :pt_in, :pt_out, :Tt_in, :Tt_out, :Mach_in, :Mach_out, :valid, :stall, :choke],
        )
    finally
        io === stdout || close(io)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
