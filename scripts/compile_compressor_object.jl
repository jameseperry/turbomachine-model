#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor

function _parse_transform(raw::AbstractString)
    key = lowercase(strip(raw))
    if key in ("spec_to_analytic_map", "spec_to_map", "spec_to_analytic")
        return :spec_to_analytic_map
    elseif key in ("design_to_spec", "design_to_compressor_spec")
        return :design_to_spec
    end
    error(
        "unsupported transform=$(raw) (expected spec_to_analytic_map|design_to_spec)",
    )
end

function _default_input_group(transform::Symbol)
    transform == :spec_to_analytic_map && return "compressor_spec"
    transform == :design_to_spec && return "compressor_design"
    error("unsupported transform=$(transform)")
end

function _default_output_group(transform::Symbol)
    transform == :spec_to_analytic_map && return "compressor_analytic_map"
    transform == :design_to_spec && return "compressor_spec"
    error("unsupported transform=$(transform)")
end

function _parse_kind(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:axial, :centrifugal) || error("kind must be axial or centrifugal")
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
        prog="compile_compressor_object.jl",
        description="Compile compressor TOML objects between design/spec/analytic-map forms.",
    )

    @add_arg_table! settings begin
        "transform"
            help = "transform: spec_to_analytic_map or design_to_spec"
            required = true
        "input_path"
            help = "input TOML path"
            required = true
        "output_path"
            help = "output TOML path"
            required = true
        "--input-group"
            help = "input TOML group/table (default depends on transform)"
            arg_type = String
        "--output-group"
            help = "output TOML group/table (default depends on transform)"
            arg_type = String
        "--kind"
            help = "compressor kind for spec_to_analytic_map: axial or centrifugal"
            arg_type = String
            default = "axial"
    end

    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    transform = _parse_transform(parsed["transform"])
    input_path = parsed["input_path"]
    output_path = parsed["output_path"]

    input_group_opt = _parsed_opt(parsed, "input_group", "input-group")
    output_group_opt = _parsed_opt(parsed, "output_group", "output-group")
    input_group = isnothing(input_group_opt) ? _default_input_group(transform) : input_group_opt
    output_group = isnothing(output_group_opt) ? _default_output_group(transform) : output_group_opt

    if transform == :spec_to_analytic_map
        spec = Compressor.read_toml(Compressor.CompressorSpec, input_path; group=input_group)
        kind = _parse_kind(parsed["kind"])
        analytic_map = Compressor.compile_compressor_map(spec; kind=kind)
        Compressor.write_toml(analytic_map, output_path; group=output_group)
        println(
            "Compiled compressor spec -> analytic map: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), kind=$(kind)",
        )
        return
    end

    if transform == :design_to_spec
        design = Compressor.read_toml(Compressor.CompressorDesign, input_path; group=input_group)
        spec = Compressor.compile_compressor_spec(design)
        Compressor.write_toml(spec, output_path; group=output_group)
        println(
            "Compiled compressor design -> spec: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group)",
        )
        return
    end

    error("unsupported transform=$(transform)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
