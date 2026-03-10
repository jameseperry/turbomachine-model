#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const TM = TurboMachineModel.Physics.Turbomachine
const AxialMachine = TM.AxialMachine

function _parse_interpolation(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:bilinear, :bicubic) ||
        error("unsupported interpolation=$(raw) (expected bilinear|bicubic)")
    return key
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _build_parser()
    settings = ArgParseSettings(
        prog="generate_performance_map_from_axial_machine_model.jl",
        description="Generate a common tabulated performance map from an axial-machine model TOML.",
    )

    @add_arg_table! settings begin
        "input_path"
            help = "input axial-machine model TOML path"
            required = true
        "output_path"
            help = "output common performance map TOML path"
            required = true
        "--input-group"
            help = "input TOML group/table name"
            arg_type = String
            default = "axial_model"
        "--output-group"
            help = "output TOML group/table name"
            arg_type = String
            default = "performance_map"
        "--interpolation"
            help = "tabulated interpolation type: bilinear or bicubic"
            arg_type = String
            default = "bilinear"
        "--n-speed"
            help = "number of corrected-speed grid points for tabulation"
            arg_type = Int
            default = 31
        "--n-flow"
            help = "number of corrected-flow grid points for tabulation"
            arg_type = Int
            default = 41
        "--boundary-resolution"
            help = "phi probe count per speed line for feasible-flow boundary detection"
            arg_type = Int
            default = 401
        "--Tt-in-ref"
            help = "reference inlet total temperature [K] for tabulation"
            arg_type = Float64
            default = 288.15
        "--Pt-in-ref"
            help = "reference inlet total pressure [Pa] for tabulation"
            arg_type = Float64
            default = 101_325.0
        "--Tt-ref"
            help = "internal corrected-speed / corrected-flow reference Tt [K]"
            arg_type = Float64
        "--Pt-ref"
            help = "internal corrected-speed / corrected-flow reference Pt [Pa]"
            arg_type = Float64
        "--omega-ref-for-phi"
            help = "reference physical omega [rad/s] used to derive the default corrected-flow grid"
            arg_type = Float64
        "--nu-theta-inlet"
            help = "inlet swirl velocity ratio for streamtube sampling"
            arg_type = Float64
            default = 0.0
        "--prefer-root"
            help = "streamtube root preference when multiple local roots exist: low or high"
            arg_type = String
            default = "low"
    end

    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    input_path = parsed["input_path"]
    output_path = parsed["output_path"]
    input_group = something(_parsed_opt(parsed, "input_group", "input-group"), "axial_model")
    output_group = something(_parsed_opt(parsed, "output_group", "output-group"), "performance_map")
    interpolation = _parse_interpolation(parsed["interpolation"])
    prefer_root = Symbol(lowercase(parsed["prefer-root"]))
    prefer_root in (:low, :high) || error("prefer-root must be low or high")

    Tt_in_ref = parsed["Tt-in-ref"]
    Pt_in_ref = parsed["Pt-in-ref"]
    Tt_ref_opt = _parsed_opt(parsed, "Tt_ref", "Tt-ref")
    Pt_ref_opt = _parsed_opt(parsed, "Pt_ref", "Pt-ref")
    Tt_ref = isnothing(Tt_ref_opt) ? Tt_in_ref : Tt_ref_opt
    Pt_ref = isnothing(Pt_ref_opt) ? Pt_in_ref : Pt_ref_opt
    omega_ref_for_phi = _parsed_opt(parsed, "omega_ref_for_phi", "omega-ref-for-phi")

    meanline_model = AxialMachine.read_toml(
        TM.AxialModel,
        input_path;
        group=input_group,
    )

    map = TM.tabulate_axial_machine_model(
        meanline_model;
        n_speed=parsed["n-speed"],
        n_flow=parsed["n-flow"],
        Tt_in_ref=Tt_in_ref,
        Pt_in_ref=Pt_in_ref,
        Tt_ref=Tt_ref,
        Pt_ref=Pt_ref,
        omega_ref_for_phi=omega_ref_for_phi,
        interpolation=interpolation,
        boundary_resolution=parsed["boundary-resolution"],
        nu_theta_inlet=parsed["nu-theta-inlet"],
        prefer_root=prefer_root,
    )
    TM.write_toml(map, output_path; group=output_group)
    n_speed = parsed["n-speed"]
    n_flow = parsed["n-flow"]
    println(
        "Generated common performance map from axial-machine model: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation), n_speed=$(n_speed), n_flow=$(n_flow), Tt_ref=$(Tt_ref), Pt_ref=$(Pt_ref), prefer_root=$(prefer_root)",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
