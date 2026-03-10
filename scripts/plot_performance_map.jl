#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots
using Plots.PlotMeasures: mm

const TM = TurboMachineModel.Physics.Turbomachine
const Fluids = TurboMachineModel.Physics.Fluids
const U = TurboMachineModel.Utility

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _load_map(path::AbstractString; group::AbstractString="performance_map")
    return TM.read_toml(TM.TabulatedPerformanceMap, path; group=group)
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

function _resample_mdot_over_omega_pr(sampled)
    finite_pr = filter(isfinite, vec(sampled.pr_eval))
    isempty(finite_pr) && error("sampled map contains no finite pressure-ratio values")
    pr_eval = collect(range(minimum(finite_pr), maximum(finite_pr), length=length(sampled.mdot_eval)))
    mdot_over_pr = fill(NaN, length(sampled.omega_eval), length(pr_eval))

    for i in eachindex(sampled.omega_eval)
        prs = Float64[]
        mdots = Float64[]
        for j in eachindex(sampled.mdot_eval)
            pr = sampled.pr_eval[i, j]
            if isfinite(pr)
                push!(prs, pr)
                push!(mdots, sampled.mdot_eval[j])
            end
        end
        length(prs) >= 2 || continue
        mdot_over_pr[i, :] .= U.resample_linear(
            prs,
            mdots,
            pr_eval;
            extrapolation=:nan,
            sort_inputs=true,
            deduplicate=:mean,
        )
    end

    return (pr_eval=pr_eval, mdot_over_pr=mdot_over_pr)
end

function plot_performance_map(
    map::TM.TabulatedPerformanceMap;
    title_prefix::String="Performance Map",
    Tt_in::Real=288.15,
    Pt_in::Real=101_325.0,
    resolution::Int=160,
)
    sampled = _sample_map(map, Tt_in, Pt_in, resolution)
    inverse_sampled = _resample_mdot_over_omega_pr(sampled)
    b = sampled.boundary

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
    p3 = _plot_contour_panel(
        sampled.omega_eval,
        inverse_sampled.pr_eval,
        inverse_sampled.mdot_over_pr';
        xlabel="omega [rad/s]",
        ylabel="PR",
        title="$(title_prefix): mdot over omega, PR",
        colorbar_title="mdot [kg/s]",
        boundary_x_low=sampled.omega_eval,
        boundary_x_high=sampled.omega_eval,
        boundary_y_low=b.pr_low,
        boundary_y_high=b.pr_high,
    )

    p4 = _plot_contour_panel(
        sampled.omega_eval,
        sampled.mdot_eval,
        sampled.power_eval';
        xlabel="omega [rad/s]",
        ylabel="mdot [kg/s]",
        title="$(title_prefix): Power over omega, mdot",
        colorbar_title="power [kW]",
        boundary_x_low=sampled.omega_eval,
        boundary_x_high=sampled.omega_eval,
        boundary_y_low=b.mdot_low,
        boundary_y_high=b.mdot_high,
    )

    return plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 1000))
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

    vals = filter(isfinite, vec(TM.tabulated_pressure_ratio_table(map)))
    title_default = isempty(vals) ? "Performance Map" : (maximum(vals) <= 1.0 ? "Turbine Map" : "Compressor Map")
    title_prefix = something(_parsed_opt(parsed, "title_prefix", "title-prefix"), title_default)
    fig = plot_performance_map(
        map;
        title_prefix=title_prefix,
        Tt_in=something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15),
        Pt_in=something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0),
        resolution=parsed["resolution"],
    )
    savefig(fig, parsed["output"])
    println("Saved map plot to: $(parsed["output"])")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
