#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots
using Plots.PlotMeasures: mm

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor
const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine

function _infer_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".toml"
        return :toml
    end
    error("unsupported map extension $(ext) for path $(path); expected .toml")
end

_default_groups(::Val{:compressor}) = ("compressor_map", "compressor_analytic_map")
_default_groups(::Val{:turbine}) = ("turbine_map",)

function _read_map(
    ::Val{:compressor},
    path::AbstractString,
    group::AbstractString,
)
    try
        return Compressor.read_toml(Compressor.TabulatedCompressorPerformanceMap, path; group=group)
    catch
    end
    return Compressor.read_toml(Compressor.AnalyticCompressorPerformanceMap, path; group=group)
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

    if kind != :auto
        groups = isnothing(group) ? _default_groups(Val(kind)) : (group,)
        for g in groups
            try
                return _read_map(Val(kind), path, g), kind
            catch
            end
        end
        error("failed to load $(kind) map from $(path); try --map-group")
    end

    for k in (:compressor, :turbine)
        groups = isnothing(group) ? _default_groups(Val(k)) : (group,)
        for g in groups
            try
                return _read_map(Val(k), path, g), k
            catch
            end
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
    if !isempty(x_points) && !isempty(y_points)
        scatter!(p_left, vec(repeat(x_points', length(y_points))), vec(repeat(y_points, 1, length(x_points))); ms=2, label=false)
    end

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
    if !isempty(x_points) && !isempty(y_points)
        scatter!(p_right, vec(repeat(x_points', length(y_points))), vec(repeat(y_points, 1, length(x_points))); ms=2, label=false)
    end

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

function _line_sample_indices(n::Int; n_lines::Int=10)
    n <= 0 && return Int[]
    k = min(n, n_lines)
    return unique(round.(Int, range(1, n, length=k)))
end

function _plot_compressor_four_panel(
    mdot_eval::AbstractVector{<:Real},
    omega_eval::AbstractVector{<:Real},
    pr_eval::AbstractMatrix{<:Real},
    eta_eval::AbstractMatrix{<:Real},
    mdot_points::AbstractVector{<:Real},
    omega_points::AbstractVector{<:Real};
    title_prefix::String,
)
    p_pr = contour(
        mdot_eval,
        omega_eval,
        pr_eval;
        xlabel="mdot_corr",
        ylabel="omega_corr",
        title="$(title_prefix): Pressure Ratio",
        colorbar_title="PR",
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    p_eta = contour(
        mdot_eval,
        omega_eval,
        eta_eval;
        xlabel="mdot_corr",
        ylabel="omega_corr",
        title="$(title_prefix): Isentropic Efficiency",
        colorbar_title="eta",
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    if !isempty(mdot_points) && !isempty(omega_points)
        xpts = vec(repeat(mdot_points', length(omega_points)))
        ypts = vec(repeat(omega_points, 1, length(mdot_points)))
        scatter!(p_pr, xpts, ypts; ms=2, label=false)
        scatter!(p_eta, xpts, ypts; ms=2, label=false)
    end

    p_pr_vs_omega = plot(
        xlabel="omega_corr",
        ylabel="PR",
        title="$(title_prefix): PR vs omega (mdot contours)",
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    p_eta_vs_omega = plot(
        xlabel="omega_corr",
        ylabel="eta",
        title="$(title_prefix): eta vs omega (mdot contours)",
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )

    for j in _line_sample_indices(length(mdot_eval))
        mdot = mdot_eval[j]
        label = "mdot=$(round(mdot; sigdigits=4))"
        plot!(p_pr_vs_omega, omega_eval, pr_eval[:, j]; lw=2, label=label)
        plot!(p_eta_vs_omega, omega_eval, eta_eval[:, j]; lw=2, label=label)
    end

    return plot(
        p_pr,
        p_eta,
        p_pr_vs_omega,
        p_eta_vs_omega;
        layout=(2, 2),
        size=(1400, 1100),
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
end

function _plot_turbine_four_panel(
    pr_eval::AbstractVector{<:Real},
    omega_eval::AbstractVector{<:Real},
    mdot_eval::AbstractMatrix{<:Real},
    eta_eval::AbstractMatrix{<:Real},
    pr_points::AbstractVector{<:Real},
    omega_points::AbstractVector{<:Real};
    title_prefix::String,
)
    p_mdot = contour(
        pr_eval,
        omega_eval,
        mdot_eval;
        xlabel="PR_turb = Pt_in/Pt_out",
        ylabel="omega_corr",
        title="$(title_prefix): Corrected Mass Flow",
        colorbar_title="mdot_corr",
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    p_eta = contour(
        pr_eval,
        omega_eval,
        eta_eval;
        xlabel="PR_turb = Pt_in/Pt_out",
        ylabel="omega_corr",
        title="$(title_prefix): Isentropic Efficiency",
        colorbar_title="eta",
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    if !isempty(pr_points) && !isempty(omega_points)
        xpts = vec(repeat(pr_points', length(omega_points)))
        ypts = vec(repeat(omega_points, 1, length(pr_points)))
        scatter!(p_mdot, xpts, ypts; ms=2, label=false)
        scatter!(p_eta, xpts, ypts; ms=2, label=false)
    end

    p_mdot_vs_omega = plot(
        xlabel="omega_corr",
        ylabel="mdot_corr",
        title="$(title_prefix): mdot vs omega (PR contours)",
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    p_eta_vs_omega = plot(
        xlabel="omega_corr",
        ylabel="eta",
        title="$(title_prefix): eta vs omega (PR contours)",
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )

    for j in _line_sample_indices(length(pr_eval))
        pr = pr_eval[j]
        label = "PR=$(round(pr; sigdigits=4))"
        plot!(p_mdot_vs_omega, omega_eval, mdot_eval[:, j]; lw=2, label=label)
        plot!(p_eta_vs_omega, omega_eval, eta_eval[:, j]; lw=2, label=label)
    end

    return plot(
        p_mdot,
        p_eta,
        p_mdot_vs_omega,
        p_eta_vs_omega;
        layout=(2, 2),
        size=(1400, 1100),
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
end

_grid_points(::Compressor.AbstractCompressorPerformanceMap) = (Float64[], Float64[])
_grid_points(map::Compressor.TabulatedCompressorPerformanceMap) = (
    collect(Compressor.mdot_corr_grid(map)),
    collect(Compressor.omega_corr_grid(map)),
)

_grid_points(::Turbine.AbstractTurbinePerformanceMap) = (Float64[], Float64[])
_grid_points(map::Turbine.TabulatedTurbinePerformanceMap) = (
    collect(Turbine.pr_turb_grid(map)),
    collect(Turbine.omega_corr_grid(map)),
)

function _sample_compressor_map(
    map::Compressor.AbstractCompressorPerformanceMap,
    resolution::Int,
)
    resolution >= 4 || error("resolution must be >= 4")
    domain = Compressor.performance_map_domain(map)
    omega_min, omega_max = domain.omega_corr
    mdot_min, mdot_max = domain.mdot_corr
    omega_eval = collect(range(omega_min, omega_max, length=resolution))
    mdot_eval = collect(range(mdot_min, mdot_max, length=resolution))
    pr_eval = Matrix{Float64}(undef, length(omega_eval), length(mdot_eval))
    eta_eval = Matrix{Float64}(undef, length(omega_eval), length(mdot_eval))
    flow_range = domain.mdot_corr_flow_range
    for (i, omega) in pairs(omega_eval)
        m_surge = flow_range.surge(omega)
        m_choke = flow_range.choke(omega)
        m_lo = min(m_surge, m_choke)
        m_hi = max(m_surge, m_choke)
        for (j, mdot) in pairs(mdot_eval)
            if mdot < m_lo || mdot > m_hi
                pr_eval[i, j] = NaN
                eta_eval[i, j] = NaN
            else
                vals = Compressor.compressor_performance_map(map, omega, mdot)
                pr_eval[i, j] = vals.PR
                eta_eval[i, j] = vals.eta
            end
        end
    end
    return mdot_eval, omega_eval, pr_eval, eta_eval
end

function plot_performance_map(
    map::Compressor.AbstractCompressorPerformanceMap;
    title_prefix::String="Compressor Map",
    resolution::Int=160,
)
    mdot_eval, omega_eval, pr_eval, eta_eval = _sample_compressor_map(map, resolution)
    mdot_points, omega_points = _grid_points(map)
    return _plot_compressor_four_panel(
        mdot_eval,
        omega_eval,
        pr_eval,
        eta_eval,
        mdot_points,
        omega_points;
        title_prefix=title_prefix,
    )
end

function plot_performance_map(
    map::Turbine.AbstractTurbinePerformanceMap;
    title_prefix::String="Turbine Map",
    resolution::Int=160,
)
    resolution >= 4 || error("resolution must be >= 4")
    domain = Turbine.performance_map_domain(map)
    omega_min, omega_max = domain.omega_corr
    pr_min, pr_max = domain.pr_turb
    pr_eval = collect(range(pr_min, pr_max, length=resolution))
    omega_eval = collect(range(omega_min, omega_max, length=resolution))
    mdot_eval = [Turbine.turbine_performance_map(map, ω, pr).mdot_corr for ω in omega_eval, pr in pr_eval]
    eta_eval = [Turbine.turbine_performance_map(map, ω, pr).eta for ω in omega_eval, pr in pr_eval]
    pr_points, omega_points = _grid_points(map)

    return _plot_turbine_four_panel(
        pr_eval,
        omega_eval,
        mdot_eval,
        eta_eval,
        pr_points,
        omega_points;
        title_prefix=title_prefix,
    )
end

function _build_parser()
    settings = ArgParseSettings(
        prog="plot_performance_map.jl",
        description="Plot compressor or turbine performance map contours.",
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
        resolution=parsed["resolution"],
    )
    savefig(fig, output_path)
    println("Saved map plot to: $(output_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
