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

function _parse_pr_grid_mode(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:union, :overlap) ||
        error("unsupported pr-grid-mode=$(raw) (expected union|overlap)")
    return key
end

function _parse_sampling_strategy(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:roots, :scan) ||
        error("unsupported sampling-strategy=$(raw) (expected roots|scan)")
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
        "--pr-grid-mode"
            help = "PR grid construction mode: union or overlap"
            arg_type = String
            default = "union"
        "--sampling-strategy"
            help = "tabulation sampling strategy: roots or scan"
            arg_type = String
            default = "roots"
        "--pr-merge-rtol"
            help = "relative tolerance for collapsing multi-valued PR scan points per speed line"
            arg_type = Float64
            default = 1e-4
        "--smooth-window"
            help = "odd moving-average window along PR axis (0/1 disables smoothing)"
            arg_type = Int
            default = 0
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
    pr_grid_mode = _parse_pr_grid_mode(parsed["pr-grid-mode"])
    sampling_strategy = _parse_sampling_strategy(parsed["sampling-strategy"])
    pr_merge_rtol = parsed["pr-merge-rtol"]
    smooth_window = parsed["smooth-window"]
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
        pr_grid_mode=pr_grid_mode,
        sampling_strategy=sampling_strategy,
        pr_merge_rtol=pr_merge_rtol,
        smooth_window=smooth_window,
    )
    Turbine.write_toml(map, output_path; group=output_group)
    println(
        "Generated turbine map from meanline model: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation), branch=$(branch), pr_grid_mode=$(pr_grid_mode), sampling_strategy=$(sampling_strategy), pr_merge_rtol=$(pr_merge_rtol), smooth_window=$(smooth_window), n_speed=$(n_speed), n_pr=$(n_pr), Tt_ref=$(Tt_ref), Pt_ref=$(Pt_ref)",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
