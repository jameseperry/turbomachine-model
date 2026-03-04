#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor

function _parse_transform(raw::AbstractString)
    key = lowercase(strip(raw))
    if key in ("dim_to_nd", "dimensional_to_nondimensional", "d2nd")
        return :dim_to_nd
    elseif key in ("nd_to_dim", "nondimensional_to_dimensional", "nd2d")
        return :nd_to_dim
    end
    error("unsupported transform=$(raw) (expected dim_to_nd|nd_to_dim)")
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _build_parser()
    settings = ArgParseSettings(
        prog="convert_compressor_map_coordinates.jl",
        description="Convert compressor performance maps between dimensional tabulated and non-dimensional tabulated TOML formats.",
    )

    @add_arg_table! settings begin
        "transform"
            help = "conversion direction: dim_to_nd or nd_to_dim"
            required = true
        "input_path"
            help = "input TOML path"
            required = true
        "output_path"
            help = "output TOML path"
            required = true
        "--input-group"
            help = "input TOML group/table name"
            arg_type = String
            default = "compressor_map"
        "--output-group"
            help = "output TOML group/table name"
            arg_type = String
            default = "compressor_map"
        "--interpolation"
            help = "output interpolation for resampled map: bilinear or bicubic"
            arg_type = String
            default = "bilinear"
        "--gamma"
            help = "ratio of specific heats for dim_to_nd conversion"
            arg_type = Float64
        "--gas-constant"
            help = "specific gas constant [J/(kg*K)] for dim_to_nd conversion"
            arg_type = Float64
        "--tip-radius-inlet"
            help = "inlet tip radius [m] for dim_to_nd conversion"
            arg_type = Float64
        "--mean-radius-inlet"
            help = "inlet mean radius [m] for dim_to_nd conversion"
            arg_type = Float64
        "--inlet-area"
            help = "inlet annulus area [m^2] for dim_to_nd conversion"
            arg_type = Float64
        "--Tt-in-ref"
            help = "reference inlet total temperature [K] used for conversion"
            arg_type = Float64
        "--Pt-in-ref"
            help = "reference inlet total pressure [Pa] used for conversion"
            arg_type = Float64
        "--Tt-ref"
            help = "output corrected-flow reference total temperature [K] for nd_to_dim conversion"
            arg_type = Float64
        "--Pt-ref"
            help = "output corrected-flow reference total pressure [Pa] for nd_to_dim conversion"
            arg_type = Float64
    end

    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    transform = _parse_transform(parsed["transform"])
    input_path = parsed["input_path"]
    output_path = parsed["output_path"]

    input_group = something(_parsed_opt(parsed, "input_group", "input-group"), "compressor_map")
    output_group = something(_parsed_opt(parsed, "output_group", "output-group"), "compressor_map")
    interpolation_raw = something(_parsed_opt(parsed, "interpolation", "interpolation"), "bilinear")
    interpolation = Symbol(lowercase(interpolation_raw))
    interpolation in (:bilinear, :bicubic) ||
        error("unsupported interpolation=$(interpolation_raw) (expected bilinear|bicubic)")

    if transform == :dim_to_nd
        gamma = _parsed_opt(parsed, "gamma", "gamma")
        gas_constant = _parsed_opt(parsed, "gas_constant", "gas-constant")
        tip_radius_inlet = _parsed_opt(parsed, "tip_radius_inlet", "tip-radius-inlet")
        mean_radius_inlet = _parsed_opt(parsed, "mean_radius_inlet", "mean-radius-inlet")
        inlet_area = _parsed_opt(parsed, "inlet_area", "inlet-area")

        isnothing(gamma) && error("--gamma is required for dim_to_nd")
        isnothing(gas_constant) && error("--gas-constant is required for dim_to_nd")
        isnothing(tip_radius_inlet) && error("--tip-radius-inlet is required for dim_to_nd")
        isnothing(mean_radius_inlet) && error("--mean-radius-inlet is required for dim_to_nd")
        isnothing(inlet_area) && error("--inlet-area is required for dim_to_nd")

        dim_map = Compressor.read_toml(Compressor.TabulatedCompressorPerformanceMap, input_path; group=input_group)

        Tt_in_ref = _parsed_opt(parsed, "Tt_in_ref", "Tt-in-ref")
        Pt_in_ref = _parsed_opt(parsed, "Pt_in_ref", "Pt-in-ref")
        nd_map = Compressor.to_nondimensional_tabulated_compressor_map(
            dim_map;
            gamma=gamma,
            gas_constant=gas_constant,
            tip_radius_inlet=tip_radius_inlet,
            mean_radius_inlet=mean_radius_inlet,
            inlet_area=inlet_area,
            Tt_in_ref=isnothing(Tt_in_ref) ? dim_map.Tt_ref : Tt_in_ref,
            Pt_in_ref=isnothing(Pt_in_ref) ? dim_map.Pt_ref : Pt_in_ref,
            interpolation=interpolation,
        )
        Compressor.write_toml(nd_map, output_path; group=output_group)
        println(
            "Converted compressor map dim->nd: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation)",
        )
        return
    end

    if transform == :nd_to_dim
        nd_map = Compressor.read_toml(
            Compressor.NonDimensionalTabulatedCompressorPerformanceMap,
            input_path;
            group=input_group,
        )

        Tt_in_ref = _parsed_opt(parsed, "Tt_in_ref", "Tt-in-ref")
        Pt_in_ref = _parsed_opt(parsed, "Pt_in_ref", "Pt-in-ref")
        Tt_ref = _parsed_opt(parsed, "Tt_ref", "Tt-ref")
        Pt_ref = _parsed_opt(parsed, "Pt_ref", "Pt-ref")
        isnothing(Tt_in_ref) && error("--Tt-in-ref is required for nd_to_dim")
        isnothing(Pt_in_ref) && error("--Pt-in-ref is required for nd_to_dim")

        dim_map = Compressor.to_tabulated_compressor_map(
            nd_map;
            Tt_in_ref=Tt_in_ref,
            Pt_in_ref=Pt_in_ref,
            Tt_ref=isnothing(Tt_ref) ? Tt_in_ref : Tt_ref,
            Pt_ref=isnothing(Pt_ref) ? Pt_in_ref : Pt_ref,
            interpolation=interpolation,
        )
        Compressor.write_toml(dim_map, output_path; group=output_group)
        println(
            "Converted compressor map nd->dim: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation)",
        )
        return
    end

    error("unsupported transform=$(transform)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
