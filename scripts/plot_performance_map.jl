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

_default_groups(::Val{:compressor}) = ("compressor_map",)
_default_groups(::Val{:turbine}) = ("turbine_map",)

function _read_map(
    ::Val{:compressor},
    path::AbstractString,
    group::AbstractString,
)
    return Compressor.read_performance_map_toml(path; group=group)
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

function _finite_range(values)
    finite = filter(isfinite, vec(values))
    isempty(finite) && return (NaN, NaN)
    return (minimum(finite), maximum(finite))
end

function _unique_sorted_pairs(xvals::AbstractVector{<:Real}, yvals::AbstractVector{<:Real})
    n = length(xvals)
    n == length(yvals) || error("xvals/yvals length mismatch")
    n == 0 && return (Float64[], Float64[])

    p = sortperm(xvals)
    xs = Float64[xvals[i] for i in p]
    ys = Float64[yvals[i] for i in p]

    x_out = Float64[]
    y_out = Float64[]
    i = 1
    while i <= n
        x = xs[i]
        j = i
        y_sum = 0.0
        c = 0
        while j <= n && xs[j] == x
            y_sum += ys[j]
            c += 1
            j += 1
        end
        push!(x_out, x)
        push!(y_out, y_sum / c)
        i = j
    end
    return (x_out, y_out)
end

function _interp_on_sorted_xy(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, x::Real)
    n = length(xs)
    n == length(ys) || error("xs/ys length mismatch")
    n >= 2 || return NaN
    x < xs[1] && return NaN
    x > xs[end] && return NaN
    i_hi = searchsortedfirst(xs, x)
    i_hi <= 1 && return ys[1]
    i_hi > n && return ys[end]
    i_lo = i_hi - 1
    xl = xs[i_lo]
    xr = xs[i_hi]
    yl = ys[i_lo]
    yr = ys[i_hi]
    xr == xl && return yl
    t = (x - xl) / (xr - xl)
    return (1 - t) * yl + t * yr
end

"""
Row-wise invert a sampled map `target_matrix(omega_i, src_j)` to build
`src(omega_i, target_k)` on a regular `(omega, target)` grid.
"""
function _rowwise_inverse_grid(
    omega_eval::AbstractVector{<:Real},
    src_eval::AbstractVector{<:Real},
    target_matrix::AbstractMatrix{<:Real},
    target_axis::AbstractVector{<:Real},
)
    size(target_matrix, 1) == length(omega_eval) || error("target_matrix row mismatch")
    size(target_matrix, 2) == length(src_eval) || error("target_matrix column mismatch")

    src_on_target = fill(NaN, length(omega_eval), length(target_axis))
    for i in eachindex(omega_eval)
        tgt_row = vec(target_matrix[i, :])
        src_row = Float64.(src_eval)
        mask = [isfinite(t) && isfinite(s) for (t, s) in zip(tgt_row, src_row)]
        count(mask) >= 2 || continue

        tgt_valid = Float64[tgt_row[k] for k in eachindex(tgt_row) if mask[k]]
        src_valid = Float64[src_row[k] for k in eachindex(src_row) if mask[k]]
        tgt_sorted, src_sorted = _unique_sorted_pairs(tgt_valid, src_valid)
        length(tgt_sorted) >= 2 || continue

        for (k, tgt) in pairs(target_axis)
            src_on_target[i, k] = _interp_on_sorted_xy(tgt_sorted, src_sorted, tgt)
        end
    end
    return src_on_target
end

function _plot_compressor_six_panel(
    mdot_eval::AbstractVector{<:Real},
    omega_eval::AbstractVector{<:Real},
    pr_eval::AbstractMatrix{<:Real},
    eta_eval::AbstractMatrix{<:Real},
    mdot_points::AbstractVector{<:Real},
    omega_points::AbstractVector{<:Real},
    boundary::NamedTuple;
    title_prefix::String,
)
    p_pr = contour(
        mdot_eval,
        omega_eval,
        pr_eval;
        xlabel="mdot (kg/s)",
        ylabel="omega (rad/s)",
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
        xlabel="mdot (kg/s)",
        ylabel="omega (rad/s)",
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
    plot!(p_pr, boundary.mdot_surge, omega_eval; lw=2, ls=:dash, color=:black, label="surge")
    plot!(p_pr, boundary.mdot_choke, omega_eval; lw=2, ls=:dot, color=:black, label="choke")
    plot!(p_eta, boundary.mdot_surge, omega_eval; lw=2, ls=:dash, color=:black, label="surge")
    plot!(p_eta, boundary.mdot_choke, omega_eval; lw=2, ls=:dot, color=:black, label="choke")

    pr_lo, pr_hi = _finite_range(pr_eval)
    eta_lo, eta_hi = _finite_range(eta_eval)
    pr_axis =
        (isfinite(pr_lo) && isfinite(pr_hi) && (pr_hi > pr_lo)) ?
        collect(range(pr_lo, pr_hi, length=120)) : Float64[]
    eta_axis =
        (isfinite(eta_lo) && isfinite(eta_hi) && (eta_hi > eta_lo)) ?
        collect(range(eta_lo, eta_hi, length=120)) : Float64[]

    mdot_on_pr = isempty(pr_axis) ? fill(NaN, length(omega_eval), 0) :
                 _rowwise_inverse_grid(omega_eval, mdot_eval, pr_eval, pr_axis)
    omega_on_pr = isempty(pr_axis) ? fill(NaN, length(mdot_eval), 0) :
                  _rowwise_inverse_grid(mdot_eval, omega_eval, permutedims(pr_eval), pr_axis)

    p_pr_vs_omega =
        isempty(pr_axis) ? plot(
            xlabel="omega (rad/s)",
            ylabel="PR",
            title="$(title_prefix): PR vs omega (mdot contours)",
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        ) : contour(
            omega_eval,
            pr_axis,
            permutedims(mdot_on_pr);
            xlabel="omega (rad/s)",
            ylabel="PR",
            title="$(title_prefix): PR vs omega (mdot contours)",
            colorbar_title="mdot (kg/s)",
            linewidth=2,
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        )
    plot!(p_pr_vs_omega, omega_eval, boundary.pr_surge; lw=2, ls=:dash, color=:black, label="surge")
    plot!(p_pr_vs_omega, omega_eval, boundary.pr_choke; lw=2, ls=:dot, color=:black, label="choke")

    mdot_on_eta = isempty(eta_axis) ? fill(NaN, length(omega_eval), 0) :
                  _rowwise_inverse_grid(omega_eval, mdot_eval, eta_eval, eta_axis)
    p_eta_vs_omega =
        isempty(eta_axis) ? plot(
            xlabel="omega (rad/s)",
            ylabel="eta",
            title="$(title_prefix): eta vs omega (mdot contours)",
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        ) : contour(
            omega_eval,
            eta_axis,
            permutedims(mdot_on_eta);
            xlabel="omega (rad/s)",
            ylabel="eta",
            title="$(title_prefix): eta vs omega (mdot contours)",
            colorbar_title="mdot (kg/s)",
            linewidth=2,
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        )
    plot!(p_eta_vs_omega, omega_eval, boundary.eta_surge; lw=2, ls=:dash, color=:black, label="surge")
    plot!(p_eta_vs_omega, omega_eval, boundary.eta_choke; lw=2, ls=:dot, color=:black, label="choke")

    p_mdot_vs_pr =
        isempty(pr_axis) ? plot(
            xlabel="PR",
            ylabel="mdot (kg/s)",
            title="$(title_prefix): mdot vs PR (omega contours)",
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        ) : contour(
            pr_axis,
            mdot_eval,
            omega_on_pr;
            xlabel="PR",
            ylabel="mdot (kg/s)",
            title="$(title_prefix): mdot vs PR (omega contours)",
            colorbar_title="omega (rad/s)",
            linewidth=2,
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        )
    plot!(p_mdot_vs_pr, boundary.pr_surge, boundary.mdot_surge; lw=2, ls=:dash, color=:black, label="surge")
    plot!(p_mdot_vs_pr, boundary.pr_choke, boundary.mdot_choke; lw=2, ls=:dot, color=:black, label="choke")

    omega_on_eta = isempty(eta_axis) ? fill(NaN, length(mdot_eval), 0) :
                   _rowwise_inverse_grid(mdot_eval, omega_eval, permutedims(eta_eval), eta_axis)
    p_mdot_vs_eta =
        isempty(eta_axis) ? plot(
            xlabel="eta",
            ylabel="mdot (kg/s)",
            title="$(title_prefix): mdot vs eta (omega contours)",
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        ) : contour(
            eta_axis,
            mdot_eval,
            omega_on_eta;
            xlabel="eta",
            ylabel="mdot (kg/s)",
            title="$(title_prefix): mdot vs eta (omega contours)",
            colorbar_title="omega (rad/s)",
            linewidth=2,
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        )
    plot!(p_mdot_vs_eta, boundary.eta_surge, boundary.mdot_surge; lw=2, ls=:dash, color=:black, label="surge")
    plot!(p_mdot_vs_eta, boundary.eta_choke, boundary.mdot_choke; lw=2, ls=:dot, color=:black, label="choke")

    return plot(
        p_pr,
        p_eta,
        p_pr_vs_omega,
        p_eta_vs_omega,
        p_mdot_vs_pr,
        p_mdot_vs_eta;
        layout=(3, 2),
        size=(1500, 1600),
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

    mdot_lo, mdot_hi = _finite_range(mdot_eval)
    eta_lo, eta_hi = _finite_range(eta_eval)
    mdot_axis =
        (isfinite(mdot_lo) && isfinite(mdot_hi) && (mdot_hi > mdot_lo)) ?
        collect(range(mdot_lo, mdot_hi, length=120)) : Float64[]
    eta_axis =
        (isfinite(eta_lo) && isfinite(eta_hi) && (eta_hi > eta_lo)) ?
        collect(range(eta_lo, eta_hi, length=120)) : Float64[]

    pr_on_mdot = isempty(mdot_axis) ? fill(NaN, length(omega_eval), 0) :
                 _rowwise_inverse_grid(omega_eval, pr_eval, mdot_eval, mdot_axis)
    pr_on_eta = isempty(eta_axis) ? fill(NaN, length(omega_eval), 0) :
                _rowwise_inverse_grid(omega_eval, pr_eval, eta_eval, eta_axis)

    p_mdot_vs_omega =
        isempty(mdot_axis) ? plot(
            xlabel="omega_corr",
            ylabel="mdot_corr",
            title="$(title_prefix): mdot vs omega (PR contours)",
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        ) : contour(
            omega_eval,
            mdot_axis,
            permutedims(pr_on_mdot);
            xlabel="omega_corr",
            ylabel="mdot_corr",
            title="$(title_prefix): mdot vs omega (PR contours)",
            colorbar_title="PR_turb",
            linewidth=2,
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        )

    p_eta_vs_omega =
        isempty(eta_axis) ? plot(
            xlabel="omega_corr",
            ylabel="eta",
            title="$(title_prefix): eta vs omega (PR contours)",
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        ) : contour(
            omega_eval,
            eta_axis,
            permutedims(pr_on_eta);
            xlabel="omega_corr",
            ylabel="eta",
            title="$(title_prefix): eta vs omega (PR contours)",
            colorbar_title="PR_turb",
            linewidth=2,
            left_margin=8mm,
            right_margin=10mm,
            bottom_margin=8mm,
            top_margin=8mm,
        )

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

_grid_points(::Turbine.AbstractTurbinePerformanceMap) = (Float64[], Float64[])
_grid_points(map::Turbine.TabulatedTurbinePerformanceMap) = (
    collect(Turbine.pr_turb_grid(map)),
    collect(Turbine.omega_corr_grid(map)),
)

function _sample_compressor_map(
    map::Compressor.AbstractCompressorPerformanceMap,
    Tt_in::Real,
    Pt_in::Real,
    resolution::Int,
)
    resolution >= 4 || error("resolution must be >= 4")
    domain = Compressor.performance_map_domain(map, Tt_in, Pt_in)
    omega_min, omega_max = domain.omega
    mdot_min_raw, mdot_max_raw = domain.mdot
    mdot_min = max(0.0, min(mdot_min_raw, mdot_max_raw))
    mdot_max = max(mdot_min, max(mdot_min_raw, mdot_max_raw))
    omega_eval = collect(range(omega_min, omega_max, length=resolution))
    mdot_eval = collect(range(mdot_min, mdot_max, length=resolution))
    pr_eval = Matrix{Float64}(undef, length(omega_eval), length(mdot_eval))
    eta_eval = Matrix{Float64}(undef, length(omega_eval), length(mdot_eval))
    flow_range = domain.mdot_flow_range
    mdot_surge = Float64[]
    mdot_choke = Float64[]
    pr_surge = Float64[]
    pr_choke = Float64[]
    eta_surge = Float64[]
    eta_choke = Float64[]

    for (i, omega) in pairs(omega_eval)
        m_surge = flow_range.surge(omega)
        m_choke = flow_range.choke(omega)
        push!(mdot_surge, m_surge)
        push!(mdot_choke, m_choke)
        vals_surge = Compressor.performance_from_stagnation(map, omega, m_surge, Tt_in, Pt_in)
        vals_choke = Compressor.performance_from_stagnation(map, omega, m_choke, Tt_in, Pt_in)
        push!(pr_surge, isfinite(vals_surge.PR) ? vals_surge.PR : NaN)
        push!(pr_choke, isfinite(vals_choke.PR) ? vals_choke.PR : NaN)
        push!(
            eta_surge,
            (isfinite(vals_surge.eta) && vals_surge.eta > 0) ? vals_surge.eta : NaN,
        )
        push!(
            eta_choke,
            (isfinite(vals_choke.eta) && vals_choke.eta > 0) ? vals_choke.eta : NaN,
        )

        m_lo = min(m_surge, m_choke)
        m_hi = max(m_surge, m_choke)
        for (j, mdot) in pairs(mdot_eval)
            if mdot < m_lo || mdot > m_hi
                pr_eval[i, j] = NaN
                eta_eval[i, j] = NaN
            else
                vals = Compressor.performance_from_stagnation(map, omega, mdot, Tt_in, Pt_in)
                pr_eval[i, j] = vals.PR
                eta_eval[i, j] = (isfinite(vals.eta) && vals.eta > 0) ? vals.eta : NaN
            end
        end
    end
    return (
        mdot_eval=mdot_eval,
        omega_eval=omega_eval,
        pr_eval=pr_eval,
        eta_eval=eta_eval,
        boundary=(
            mdot_surge=mdot_surge,
            mdot_choke=mdot_choke,
            pr_surge=pr_surge,
            pr_choke=pr_choke,
            eta_surge=eta_surge,
            eta_choke=eta_choke,
        ),
    )
end

function plot_performance_map(
    map::Compressor.AbstractCompressorPerformanceMap;
    title_prefix::String="Compressor Map",
    Tt_in::Real=288.15,
    Pt_in::Real=101_325.0,
    resolution::Int=160,
)
    sampled = _sample_compressor_map(map, Tt_in, Pt_in, resolution)
    mdot_points, omega_points = _grid_points(map)
    return _plot_compressor_six_panel(
        sampled.mdot_eval,
        sampled.omega_eval,
        sampled.pr_eval,
        sampled.eta_eval,
        mdot_points,
        omega_points,
        sampled.boundary;
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
    eta_eval = [
        begin
            eta = Turbine.turbine_performance_map(map, ω, pr).eta
            (isfinite(eta) && eta > 0) ? eta : NaN
        end for ω in omega_eval, pr in pr_eval
    ]
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
        "--tt-in"
            help = "inlet total temperature (K) for compressor physical-domain plotting"
            arg_type = Float64
            default = 288.15
        "--pt-in"
            help = "inlet total pressure (Pa) for compressor physical-domain plotting"
            arg_type = Float64
            default = 101_325.0
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

    fig = if detected_kind == :compressor
        tt_in = something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15)
        pt_in = something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0)
        plot_performance_map(
            map;
            title_prefix=title_prefix,
            Tt_in=tt_in,
            Pt_in=pt_in,
            resolution=parsed["resolution"],
        )
    else
        plot_performance_map(
            map;
            title_prefix=title_prefix,
            resolution=parsed["resolution"],
        )
    end
    savefig(fig, output_path)
    println("Saved map plot to: $(output_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
