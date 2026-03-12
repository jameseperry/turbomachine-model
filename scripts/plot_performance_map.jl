#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots
using Plots.PlotMeasures: mm

const TM = TurboMachineModel.Physics.Turbomachine
const Fluids = TurboMachineModel.Physics.Fluids

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _load_map(path::AbstractString; group::AbstractString="performance_map")
    return TM.read_toml(TM.TabulatedPerformanceMap, path; group=group)
end

function _load_diagnostics(path::AbstractString; group::AbstractString="performance_map_diagnostics")
    return TM.read_toml(TM.TabulatedPerformanceMapDiagnostics, path; group=group)
end

function _sample_map(
    map::TM.TabulatedPerformanceMap,
    Tt_in::Real,
    Pt_in::Real,
    resolution::Int,
)
    resolution >= 4 || error("resolution must be >= 4")
    domain = TM.performance_map_domain(map, Tt_in, Pt_in)
    omega_min, omega_max = domain.omega
    mdot_min_raw, mdot_max_raw = domain.mdot
    mdot_min = max(0.0, min(mdot_min_raw, mdot_max_raw))
    mdot_max = max(mdot_min, max(mdot_min_raw, mdot_max_raw))

    omega_eval = collect(range(omega_min, omega_max, length=resolution))
    mdot_eval = collect(range(mdot_min, mdot_max, length=resolution))
    pr_eval = fill(NaN, length(omega_eval), length(mdot_eval))
    eta_eval = fill(NaN, length(omega_eval), length(mdot_eval))
    power_eval = fill(NaN, length(omega_eval), length(mdot_eval))

    eos = Fluids.ideal_EOS()[:air]
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)

    mdot_low = Float64[]
    mdot_high = Float64[]
    pr_low = Float64[]
    pr_high = Float64[]
    eta_low = Float64[]
    eta_high = Float64[]

    for (i, omega) in pairs(omega_eval)
        m_low = domain.mdot_flow_range.min(omega)
        m_high = domain.mdot_flow_range.max(omega)
        push!(mdot_low, m_low)
        push!(mdot_high, m_high)

        vals_low = TM.performance_from_stagnation(map, omega, m_low, Tt_in, Pt_in)
        vals_high = TM.performance_from_stagnation(map, omega, m_high, Tt_in, Pt_in)
        push!(pr_low, isfinite(vals_low.pressure_ratio) ? vals_low.pressure_ratio : NaN)
        push!(pr_high, isfinite(vals_high.pressure_ratio) ? vals_high.pressure_ratio : NaN)
        push!(eta_low, (isfinite(vals_low.eta) && vals_low.eta > 0) ? vals_low.eta : NaN)
        push!(eta_high, (isfinite(vals_high.eta) && vals_high.eta > 0) ? vals_high.eta : NaN)

        m_lo = min(m_low, m_high)
        m_hi = max(m_low, m_high)
        for (j, mdot) in pairs(mdot_eval)
            if m_lo <= mdot <= m_hi
                vals = TM.performance_from_stagnation(map, omega, mdot, Tt_in, Pt_in)
                pr_eval[i, j] = vals.pressure_ratio
                eta_eval[i, j] = (isfinite(vals.eta) && vals.eta > 0) ? vals.eta : NaN
                if isfinite(vals.pressure_ratio) && isfinite(vals.eta) && vals.eta > 0
                    pt_out = vals.pressure_ratio * Pt_in
                    s_in = Fluids.entropy(eos, Pt_in, ht_in)
                    h2s = Fluids.enthalpy_from_pressure_entropy(eos, pt_out, s_in)
                    if h2s >= ht_in
                        ht_out = ht_in + (h2s - ht_in) / vals.eta
                    else
                        ht_out = ht_in - vals.eta * (ht_in - h2s)
                    end
                    power_eval[i, j] = mdot * (ht_out - ht_in) / 1_000.0
                end
            end
        end
    end

    return (
        mdot_eval=mdot_eval,
        omega_eval=omega_eval,
        pr_eval=pr_eval,
        eta_eval=eta_eval,
        power_eval=power_eval,
        boundary=(
            mdot_low=mdot_low,
            mdot_high=mdot_high,
            pr_low=pr_low,
            pr_high=pr_high,
            eta_low=eta_low,
            eta_high=eta_high,
        ),
    )
end

function _plot_contour_panel(
    x,
    y,
    z;
    xlabel::String,
    ylabel::String,
    title::String,
    colorbar_title::String,
    boundary_x_low::AbstractVector{<:Real}=Float64[],
    boundary_x_high::AbstractVector{<:Real}=Float64[],
    boundary_y_low::AbstractVector{<:Real}=Float64[],
    boundary_y_high::AbstractVector{<:Real}=Float64[],
)
    p = contour(
        x,
        y,
        z;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        colorbar_title=colorbar_title,
        linewidth=2,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    if !isempty(boundary_x_low)
        plot!(p, boundary_x_low, boundary_y_low; lw=2, ls=:dash, color=:black, label="low-flow")
    end
    if !isempty(boundary_x_high)
        plot!(p, boundary_x_high, boundary_y_high; lw=2, ls=:dot, color=:black, label="high-flow")
    end
    return p
end

function _plot_modes_panel(
    ;
    xlabel::String,
    ylabel::String,
    title::String,
    boundary_x_low::AbstractVector{<:Real}=Float64[],
    boundary_x_high::AbstractVector{<:Real}=Float64[],
    boundary_y_low::AbstractVector{<:Real}=Float64[],
    boundary_y_high::AbstractVector{<:Real}=Float64[],
)
    p = plot(
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        legend=:outertopright,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    if !isempty(boundary_x_low)
        plot!(p, boundary_x_low, boundary_y_low; lw=2, ls=:dash, color=:black, label="low-flow")
    end
    if !isempty(boundary_x_high)
        plot!(p, boundary_x_high, boundary_y_high; lw=2, ls=:dot, color=:black, label="high-flow")
    end
    return p
end

function _sample_regime(sample)
    result = sample.result
    !result.valid && return :invalid
    for row in result.row_data
        row.aero.regime == :normal || return row.aero.regime
    end
    return :normal
end

function _diagnostic_overlay_points(diagnostics::TM.TabulatedPerformanceMapDiagnostics)
    points = NamedTuple[]
    for sample in diagnostics.samples
        result = sample.result
        isfinite(result.pressure_ratio) || continue
        push!(points, (
            omega=sample.omega,
            mdot=sample.mdot_used,
            PR=result.pressure_ratio,
            regime=_sample_regime(sample),
        ))
    end
    return points
end

function _sample_first_rotor_row_incidence_deg(sample)
    rows = sample.result.row_data
    stations = sample.result.stations
    isempty(rows) && return NaN
    for row in rows
        i = row.row_index
        i + 1 <= length(stations) || continue
        delta_ht = stations[i + 1].ht_t - stations[i].ht_t
        abs(delta_ht) > 1e-8 && return row.aero.incidence * 180 / pi
    end
    return first(rows).aero.incidence * 180 / pi
end

function _diagnostic_incidence_points(diagnostics::TM.TabulatedPerformanceMapDiagnostics)
    points = NamedTuple[]
    for sample in diagnostics.samples
        result = sample.result
        incidence_deg = _sample_first_rotor_row_incidence_deg(sample)
        isfinite(result.pressure_ratio) || continue
        isfinite(incidence_deg) || continue
        push!(points, (
            omega=sample.omega,
            mdot=sample.mdot_used,
            PR=result.pressure_ratio,
            incidence_deg=incidence_deg,
        ))
    end
    return points
end

function _diagnostic_incidence_grids(diagnostics::TM.TabulatedPerformanceMapDiagnostics)
    return diagnostics
end

function _sample_incidence_over_omega_mdot(
    map::TM.TabulatedPerformanceMap,
    diagnostics::TM.TabulatedPerformanceMapDiagnostics,
    Tt_in::Real,
    Pt_in::Real,
    omega_eval::AbstractVector{<:Real},
    mdot_eval::AbstractVector{<:Real},
)
    incidence_eval = fill(NaN, length(omega_eval), length(mdot_eval))
    for (i, omega) in pairs(omega_eval), (j, mdot) in pairs(mdot_eval)
        idx = TM.nearest_grid_index(map, omega, mdot, Tt_in, Pt_in)
        incidence_eval[i, j] = _sample_first_rotor_row_incidence_deg(diagnostics.samples[idx])
    end
    return incidence_eval
end

function _resample_scalar_over_omega_pr(sampled, scalar_eval)
    finite_pr = filter(isfinite, vec(sampled.pr_eval))
    isempty(finite_pr) && error("sampled map contains no finite pressure-ratio values")
    pr_eval = collect(range(minimum(finite_pr), maximum(finite_pr), length=length(sampled.mdot_eval)))
    scalar_over_pr = fill(NaN, length(sampled.omega_eval), length(pr_eval))

    for i in eachindex(sampled.omega_eval)
        prs = Float64[]
        vals = Float64[]
        for j in eachindex(sampled.mdot_eval)
            pr = sampled.pr_eval[i, j]
            val = scalar_eval[i, j]
            if isfinite(pr) && isfinite(val)
                push!(prs, pr)
                push!(vals, val)
            end
        end
        length(prs) >= 2 || continue
        scalar_over_pr[i, :] .= TurboMachineModel.Utility.resample_linear(
            prs,
            vals,
            pr_eval;
            extrapolation=:nan,
            sort_inputs=true,
            deduplicate=:mean,
        )
    end

    return (pr_eval=pr_eval, values=scalar_over_pr)
end

function _overlay_sampling_points!(plt, points; yfield::Symbol, markersize::Real=4)
    isempty(points) && return plt
    palette = Dict(
        :normal => :forestgreen,
        :stall => :darkorange,
        :invalid => :firebrick,
    )
    labels = Dict(
        :normal => "sampling: normal",
        :stall => "sampling: stall",
        :invalid => "sampling: invalid",
    )
    for regime in (:normal, :stall, :invalid)
        regime_points = filter(p -> p.regime == regime, points)
        isempty(regime_points) && continue
        scatter!(
            plt,
            [p.omega for p in regime_points],
            [getproperty(p, yfield) for p in regime_points];
            color=palette[regime],
            markersize=markersize,
            markerstrokewidth=0,
            alpha=0.9,
            label=labels[regime],
        )
    end
    return plt
end

function _overlay_incidence_points!(plt, points; yfield::Symbol, markersize::Real=4)
    isempty(points) && return plt
    scatter!(
        plt,
        [p.omega for p in points],
        [getproperty(p, yfield) for p in points];
        marker_z=[p.incidence_deg for p in points],
        color=:viridis,
        colorbar_title="max |incidence| [deg]",
        markersize=markersize,
        markerstrokewidth=0,
        alpha=0.9,
        label="",
    )
    return plt
end

function plot_performance_map(
    map::TM.TabulatedPerformanceMap;
    title_prefix::String="Performance Map",
    Tt_in::Real=288.15,
    Pt_in::Real=101_325.0,
    resolution::Int=160,
    overlay_points=NamedTuple[],
    incidence_grids=nothing,
)
    sampled = _sample_map(map, Tt_in, Pt_in, resolution)
    b = sampled.boundary

    isempty(overlay_points) && error("plot_performance_map requires diagnostics overlay points")
    isnothing(incidence_grids) && error("plot_performance_map requires diagnostics incidence grids")
    incidence_eval = _sample_incidence_over_omega_mdot(
        map,
        incidence_grids,
        Tt_in,
        Pt_in,
        sampled.omega_eval,
        sampled.mdot_eval,
    )
    incidence_pr = _resample_scalar_over_omega_pr(sampled, incidence_eval)

    p1 = _plot_contour_panel(
        sampled.omega_eval,
        sampled.mdot_eval,
        sampled.pr_eval';
        xlabel="omega [rad/s]",
        ylabel="mdot [kg/s]",
        title="$(title_prefix): PR over omega, mdot",
        colorbar_title="PR",
        boundary_x_low=sampled.omega_eval,
        boundary_x_high=sampled.omega_eval,
        boundary_y_low=b.mdot_low,
        boundary_y_high=b.mdot_high,
    )
    _overlay_sampling_points!(p1, overlay_points; yfield=:mdot, markersize=2)
    p2 = _plot_contour_panel(
        sampled.omega_eval,
        sampled.mdot_eval,
        sampled.eta_eval';
        xlabel="omega [rad/s]",
        ylabel="mdot [kg/s]",
        title="$(title_prefix): eta over omega, mdot",
        colorbar_title="eta",
        boundary_x_low=sampled.omega_eval,
        boundary_x_high=sampled.omega_eval,
        boundary_y_low=b.mdot_low,
        boundary_y_high=b.mdot_high,
    )
    _overlay_sampling_points!(p2, overlay_points; yfield=:mdot, markersize=2)

    p3 = _plot_modes_panel(
        xlabel="omega [rad/s]",
        ylabel="mdot [kg/s]",
        title="$(title_prefix): Operating Modes over omega, mdot",
        boundary_x_low=sampled.omega_eval,
        boundary_x_high=sampled.omega_eval,
        boundary_y_low=b.mdot_low,
        boundary_y_high=b.mdot_high,
    )
    _overlay_sampling_points!(p3, overlay_points; yfield=:mdot)
    p4 = _plot_modes_panel(
        xlabel="omega [rad/s]",
        ylabel="PR",
        title="$(title_prefix): Operating Modes over omega, PR",
        boundary_x_low=sampled.omega_eval,
        boundary_x_high=sampled.omega_eval,
        boundary_y_low=b.pr_low,
        boundary_y_high=b.pr_high,
    )
    _overlay_sampling_points!(p4, overlay_points; yfield=:PR)

    p5 = contourf(
        sampled.omega_eval,
        sampled.mdot_eval,
        incidence_eval';
        xlabel="omega [rad/s]",
        ylabel="mdot [kg/s]",
        title="$(title_prefix): First-Rotor Incidence over omega, mdot",
        colorbar_title="incidence [deg]",
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    plot!(p5, sampled.omega_eval, b.mdot_low; lw=2, ls=:dash, color=:black, label="")
    plot!(p5, sampled.omega_eval, b.mdot_high; lw=2, ls=:dot, color=:black, label="")
    p6 = contourf(
        sampled.omega_eval,
        incidence_pr.pr_eval,
        incidence_pr.values';
        xlabel="omega [rad/s]",
        ylabel="PR",
        title="$(title_prefix): First-Rotor Incidence over omega, PR",
        colorbar_title="incidence [deg]",
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )
    plot!(p6, sampled.omega_eval, b.pr_low; lw=2, ls=:dash, color=:black, label="")
    plot!(p6, sampled.omega_eval, b.pr_high; lw=2, ls=:dot, color=:black, label="")

    return plot(p1, p2, p3, p4, p5, p6; layout=(3, 2), size=(1400, 1500))
end

function _build_parser()
    settings = ArgParseSettings(
        prog="plot_performance_map.jl",
        description="Plot contours for a common turbomachine performance map.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input performance map (.toml)"
            required = true
        "--map-group"
            help = "input map group/table"
            arg_type = String
            default = "performance_map"
        "--diagnostics-path"
            help = "Diagnostics TOML path on the same sampled grid as the map."
            arg_type = String
            required = true
        "--diagnostics-group"
            help = "input diagnostics group/table"
            arg_type = String
            default = "performance_map_diagnostics"
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
            help = "inlet total temperature (K)"
            arg_type = Float64
            default = 288.15
        "--pt-in"
            help = "inlet total pressure (Pa)"
            arg_type = Float64
            default = 101_325.0
    end
    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    map_group = String(something(_parsed_opt(parsed, "map_group", "map-group"), "performance_map"))
    map = _load_map(parsed["map_path"]; group=map_group)
    diagnostics_path = _parsed_opt(parsed, "diagnostics_path", "diagnostics-path")
    diagnostics_group = String(something(_parsed_opt(parsed, "diagnostics_group", "diagnostics-group"), "performance_map_diagnostics"))
    diagnostics = _load_diagnostics(String(diagnostics_path); group=diagnostics_group)
    overlay_points = _diagnostic_overlay_points(diagnostics)
    incidence_grids = _diagnostic_incidence_grids(diagnostics)

    vals = filter(isfinite, vec(TM.tabulated_pressure_ratio_table(map)))
    title_default = isempty(vals) ? "Performance Map" : (maximum(vals) <= 1.0 ? "Turbine Map" : "Compressor Map")
    title_prefix = something(_parsed_opt(parsed, "title_prefix", "title-prefix"), title_default)
    fig = plot_performance_map(
        map;
        title_prefix=title_prefix,
        Tt_in=something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15),
        Pt_in=something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0),
        resolution=parsed["resolution"],
        overlay_points=overlay_points,
        incidence_grids=incidence_grids,
    )
    savefig(fig, parsed["output"])
    println("Saved map plot to: $(parsed["output"])")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
