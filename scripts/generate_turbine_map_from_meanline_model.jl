#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine

function _parse_interpolation(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:bilinear, :bicubic) ||
        error("unsupported interpolation=$(raw) (expected bilinear|bicubic)")
    return key
end

function _parse_branch(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:low, :high) || error("unsupported branch=$(raw) (expected low|high)")
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
        prog="generate_turbine_map_from_meanline_model.jl",
        description="Generate a tabulated turbine performance map from an axial-machine meanline model TOML.",
    )

    @add_arg_table! settings begin
        "input_path"
            help = "input meanline model TOML path"
            required = true
        "output_path"
            help = "output turbine map TOML path"
            required = true
        "--input-group"
            help = "input TOML group/table name"
            arg_type = String
            default = "compressor_meanline_model"
        "--output-group"
            help = "output TOML group/table name"
            arg_type = String
            default = "turbine_map"
        "--interpolation"
            help = "tabulated interpolation type: bilinear or bicubic"
            arg_type = String
            default = "bilinear"
        "--branch"
            help = "phi-root branch selection when multiple solutions exist: low or high"
            arg_type = String
            default = "high"
        "--n-speed"
            help = "number of speed grid points for tabulation"
            arg_type = Int
            default = 31
        "--n-pr"
            help = "number of PR_turb grid points for tabulation"
            arg_type = Int
            default = 41
        "--boundary-resolution"
            help = "phi probe count per speed line for feasible PR range detection/root bracketing"
            arg_type = Int
            default = 401
        "--Tt-in-ref"
            help = "reference inlet total temperature [K]"
            arg_type = Float64
            default = 288.15
        "--Pt-in-ref"
            help = "reference inlet total pressure [Pa]"
            arg_type = Float64
            default = 101_325.0
        "--Tt-ref"
            help = "output corrected-flow reference Tt [K]"
            arg_type = Float64
        "--Pt-ref"
            help = "output corrected-flow reference Pt [Pa]"
            arg_type = Float64
    end

    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    input_path = parsed["input_path"]
    output_path = parsed["output_path"]
    input_group = something(_parsed_opt(parsed, "input_group", "input-group"), "compressor_meanline_model")
    output_group = something(_parsed_opt(parsed, "output_group", "output-group"), "turbine_map")
    interpolation = _parse_interpolation(parsed["interpolation"])
    branch = _parse_branch(parsed["branch"])
    n_speed = parsed["n-speed"]
    n_pr = parsed["n-pr"]
    boundary_resolution = parsed["boundary-resolution"]
    Tt_in_ref = parsed["Tt-in-ref"]
    Pt_in_ref = parsed["Pt-in-ref"]
    Tt_ref_opt = _parsed_opt(parsed, "Tt_ref", "Tt-ref")
    Pt_ref_opt = _parsed_opt(parsed, "Pt_ref", "Pt-ref")
    Tt_ref = isnothing(Tt_ref_opt) ? Tt_in_ref : Tt_ref_opt
    Pt_ref = isnothing(Pt_ref_opt) ? Pt_in_ref : Pt_ref_opt

    meanline_model = Turbine.read_toml(
        Turbine.TurbineMeanlineModel,
        input_path;
        group=input_group,
    )
    map = Turbine.tabulate_turbine_meanline_model(
        meanline_model;
        n_speed=n_speed,
        n_pr=n_pr,
        interpolation=interpolation,
        boundary_resolution=boundary_resolution,
        Tt_in_ref=Tt_in_ref,
        Pt_in_ref=Pt_in_ref,
        Tt_ref=Tt_ref,
        Pt_ref=Pt_ref,
        branch=branch,
    )
    Turbine.write_toml(map, output_path; group=output_group)
    println(
        "Generated turbine map from meanline model: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation), branch=$(branch), n_speed=$(n_speed), n_pr=$(n_pr), Tt_ref=$(Tt_ref), Pt_ref=$(Pt_ref)",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
