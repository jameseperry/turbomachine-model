#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots
using Plots.PlotMeasures: mm

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor
const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine
const U = TurboMachineModel.Utility

function _infer_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".toml"
        return :toml
    end
    error("unsupported map extension $(ext) for path $(path); expected .toml")
end

_default_group(::Val{:compressor}) = "compressor_map"
_default_group(::Val{:turbine}) = "turbine_map"

function _read_map(
    ::Val{:compressor},
    path::AbstractString,
    group::AbstractString,
)
    return Compressor.read_toml(Compressor.TabulatedCompressorPerformanceMap, path; group=group)
end

function _read_map(
    ::Val{:turbine},
    path::AbstractString,
    group::AbstractString,
)
    return Turbine.read_toml(Turbine.TabulatedTurbinePerformanceMap, path; group=group)
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _load_map(path::AbstractString; kind::Symbol=:auto, group::Union{Nothing,String}=nothing)
    _infer_format(path)
    kind in (:auto, :compressor, :turbine) ||
        error("unsupported kind=$(kind), expected auto|compressor|turbine")

    if kind == :compressor
        g = isnothing(group) ? _default_group(Val(:compressor)) : group
        return _read_map(Val(:compressor), path, g), :compressor
    elseif kind == :turbine
        g = isnothing(group) ? _default_group(Val(:turbine)) : group
        return _read_map(Val(:turbine), path, g), :turbine
    end

    for k in (:compressor, :turbine)
        g = isnothing(group) ? _default_group(Val(k)) : group
        try
            return _read_map(Val(k), path, g), k
        catch
        end
    end

    error("failed to auto-load map from $(path); specify --kind and/or --map-group")
end

function _plot_two_panel(
    x_eval::AbstractVector{<:Real},
    y_eval::AbstractVector{<:Real},
    z_left::AbstractMatrix{<:Real},
    z_right::AbstractMatrix{<:Real},
    x_points::AbstractVector{<:Real},
    y_points::AbstractVector{<:Real};
    x_label::String,
    y_label::String,
    left_title::String,
    right_title::String,
    left_colorbar::String,
    right_colorbar::String,
)
    p_left = contour(
        x_eval,
        y_eval,
        z_left;
        xlabel=x_label,
        ylabel=y_label,
        title=left_title,
        colorbar_title=left_colorbar,
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    scatter!(p_left, vec(repeat(x_points', length(y_points))), vec(repeat(y_points, 1, length(x_points))); ms=2, label=false)

    p_right = contour(
        x_eval,
        y_eval,
        z_right;
        xlabel=x_label,
        ylabel=y_label,
        title=right_title,
        colorbar_title=right_colorbar,
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    scatter!(p_right, vec(repeat(x_points', length(y_points))), vec(repeat(y_points, 1, length(x_points))); ms=2, label=false)

    return plot(
        p_left,
        p_right;
        layout=(1, 2),
        size=(1300, 650),
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
end

function plot_performance_map(
    map::Compressor.TabulatedCompressorPerformanceMap;
    title_prefix::String="Compressor Map",
    smooth::Bool=true,
    resolution::Int=160,
)
    mdot_grid = U.table_ygrid(map.pr_map)
    omega_grid = U.table_xgrid(map.pr_map)
    pr_grid = U.table_values(map.pr_map)
    eta_grid = U.table_values(map.eta_map)

    mdot_eval = mdot_grid
    omega_eval = omega_grid
    pr_eval = pr_grid
    eta_eval = eta_grid

    if smooth
        resolution >= 4 || error("resolution must be >= 4")
        mdot_eval = collect(range(first(mdot_grid), last(mdot_grid), length=resolution))
        omega_eval = collect(range(first(omega_grid), last(omega_grid), length=resolution))
        pr_eval = [Compressor.compressor_performance_map(map, ω, m).PR for ω in omega_eval, m in mdot_eval]
        eta_eval = [Compressor.compressor_performance_map(map, ω, m).eta for ω in omega_eval, m in mdot_eval]
    end

    return _plot_two_panel(
        mdot_eval,
        omega_eval,
        pr_eval,
        eta_eval,
        mdot_grid,
        omega_grid;
        x_label="mdot_corr",
        y_label="omega_corr",
        left_title="$(title_prefix): Pressure Ratio",
        right_title="$(title_prefix): Isentropic Efficiency",
        left_colorbar="PR",
        right_colorbar="eta",
    )
end

function plot_performance_map(
    map::Turbine.TabulatedTurbinePerformanceMap;
    title_prefix::String="Turbine Map",
    smooth::Bool=true,
    resolution::Int=160,
)
    pr_grid = U.table_ygrid(map.mdot_corr_map)
    omega_grid = U.table_xgrid(map.mdot_corr_map)
    mdot_grid = U.table_values(map.mdot_corr_map)
    eta_grid = U.table_values(map.eta_map)

    pr_eval = pr_grid
    omega_eval = omega_grid
    mdot_eval = mdot_grid
    eta_eval = eta_grid

    if smooth
        resolution >= 4 || error("resolution must be >= 4")
        pr_eval = collect(range(first(pr_grid), last(pr_grid), length=resolution))
        omega_eval = collect(range(first(omega_grid), last(omega_grid), length=resolution))
        mdot_eval = [Turbine.turbine_performance_map(map, ω, pr).mdot_corr for ω in omega_eval, pr in pr_eval]
        eta_eval = [Turbine.turbine_performance_map(map, ω, pr).eta for ω in omega_eval, pr in pr_eval]
    end

    return _plot_two_panel(
        pr_eval,
        omega_eval,
        mdot_eval,
        eta_eval,
        pr_grid,
        omega_grid;
        x_label="PR_turb = Pt_in/Pt_out",
        y_label="omega_corr",
        left_title="$(title_prefix): Corrected Mass Flow",
        right_title="$(title_prefix): Isentropic Efficiency",
        left_colorbar="mdot_corr",
        right_colorbar="eta",
    )
end

function _build_parser()
    settings = ArgParseSettings(
        prog="plot_performance_map.jl",
        description="Plot compressor or turbine tabulated performance map contours.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input performance map (.toml)"
            required = true
        "--kind"
            help = "map kind: auto, compressor, or turbine"
            arg_type = String
            default = "auto"
        "--map-group"
            help = "input map group/table"
            arg_type = String
        "--output"
            help = "output plot path"
            arg_type = String
            default = "performance_map.png"
        "--title-prefix"
            help = "plot title prefix"
            arg_type = String
        "--resolution"
            help = "grid resolution for interpolated plotting"
            arg_type = Int
            default = 160
    end
    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    kind_raw = lowercase(parsed["kind"])
    kind = Symbol(kind_raw)
    kind in (:auto, :compressor, :turbine) ||
        error("unsupported kind=$(kind_raw) (expected auto|compressor|turbine)")
    map_group = _parsed_opt(parsed, "map_group", "map-group")
    map, detected_kind = _load_map(parsed["map_path"]; kind=kind, group=map_group)

    title_default = detected_kind == :compressor ? "Compressor Map" : "Turbine Map"
    title_prefix = something(_parsed_opt(parsed, "title_prefix", "title-prefix"), title_default)
    output_path = parsed["output"]

    fig = plot_performance_map(
        map;
        title_prefix=title_prefix,
        smooth=true,
        resolution=parsed["resolution"],
    )
    savefig(fig, output_path)
    println("Saved map plot to: $(output_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
