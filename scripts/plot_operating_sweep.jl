#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots
using Plots.PlotMeasures: mm

const TM = TurboMachineModel.Physics.Turbomachine
const Fluids = TurboMachineModel.Physics.Fluids

function _parsed_opt(parsed::Dict, underscored::AbstractString, hyphenated::AbstractString; default=nothing)
    if haskey(parsed, underscored)
        return parsed[underscored]
    elseif haskey(parsed, hyphenated)
        return parsed[hyphenated]
    elseif default !== nothing
        return default
    end
    error("missing required parsed option: $(underscored) / $(hyphenated)")
end

function _load_map(path::AbstractString; group::AbstractString="performance_map")
    return TM.read_toml(TM.TabulatedPerformanceMap, path; group=group)
end

function _load_diagnostics(path::AbstractString; group::AbstractString="performance_map_diagnostics")
    return TM.read_toml(TM.TabulatedPerformanceMapDiagnostics, path; group=group)
end

function _default_omega_bounds(map::TM.TabulatedPerformanceMap, Tt_in::Real, pt_in::Real)
    domain = TM.performance_map_domain(map, Tt_in, pt_in)
    return domain.omega
end

function _default_pr_bounds(map::TM.TabulatedPerformanceMap, machine_kind::Symbol)
    vals = vec(TM.tabulated_pressure_ratio_table(map))
    finite_vals = filter(isfinite, vals)
    if machine_kind == :compressor
        lo = max(1.0, minimum(finite_vals))
        hi = max(lo, maximum(finite_vals))
        return (lo, hi)
    elseif machine_kind == :turbine
        lo = min(0.95, minimum(finite_vals))
        hi = min(1.0, maximum(finite_vals))
        return (lo, hi)
    end
    error("machine_kind must be :compressor or :turbine")
end

function _resolve_grid(min_val::Union{Nothing,Real}, max_val::Union{Nothing,Real}, n::Int, default_bounds)
    n >= 2 || error("grid resolution must be >= 2")
    lo = isnothing(min_val) ? default_bounds[1] : Float64(min_val)
    hi = isnothing(max_val) ? default_bounds[2] : Float64(max_val)
    hi > lo || error("grid max must be > grid min")
    return collect(range(lo, hi, length=n))
end

function _branch_count_matrix(point_results)
    n1, n2 = size(point_results)
    counts = Matrix{Float64}(undef, n1, n2)
    for i in 1:n1, j in 1:n2
        counts[i, j] = length(point_results[i, j].candidates)
    end
    return counts
end

function _candidate_aggregate_matrix(point_results, f::Function)
    n1, n2 = size(point_results)
    vals = fill(NaN, n1, n2)
    for i in 1:n1, j in 1:n2
        candidates = point_results[i, j].candidates
        isempty(candidates) && continue
        vals[i, j] = f(candidates)
    end
    return vals
end

function _candidate_extreme_field_matrix(point_results, field::Symbol, which::Symbol)
    which in (:low, :high) || error("which must be :low or :high")
    return _candidate_aggregate_matrix(point_results, candidates -> begin
        target = which == :low ? argmin(c.mdot for c in candidates) : argmax(c.mdot for c in candidates)
        getproperty(candidates[target], field)
    end)
end

function _single_rows(data)
    rows = NamedTuple[]
    for i in eachindex(data.omega_grid)
        for j in eachindex(data.pt_out_grid)
            push!(rows, (
                i_omega=i,
                i_pt_out=j,
                omega=data.omega_grid[i],
                pt_out=data.pt_out_grid[j],
                PR=data.PR[i, j],
                mdot=data.mdot[i, j],
                ht_out=data.ht_out[i, j],
                tau=data.tau[i, j],
                power=data.power[i, j],
                eta=data.eta[i, j],
                converged=data.converged[i, j],
                branch_coordinate=data.branch_coordinate[i, j],
                branch_index=data.branch_index[i, j],
            ))
        end
    end
    return rows
end

_all_rows(data) = data.rows

function _write_csv(path::AbstractString, rows::Vector{<:NamedTuple})
    isempty(rows) && error("cannot write CSV with no rows")
    headers = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(String.(headers), ","))
        for row in rows
            vals = [getproperty(row, h) for h in headers]
            println(io, join(map(v -> string(v), vals), ","))
        end
    end
    return path
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
            PR=result.pressure_ratio,
            regime=_sample_regime(sample),
        ))
    end
    return points
end

function _overlay_sampling_points!(plt, points)
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
            [p.PR for p in regime_points];
            color=palette[regime],
            markersize=2.5,
            markerstrokewidth=0,
            alpha=0.8,
            label=labels[regime],
        )
    end
    return plt
end

function _plot_single(data; output_path::AbstractString, overlay_points=NamedTuple[])
    pt_in = data.diagnostics[1, 1].pt_in
    pr_grid = data.pt_out_grid ./ pt_in
    branch_count = _branch_count_matrix(data.point_results)

    p1 = heatmap(data.omega_grid, pr_grid, permutedims(data.mdot); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Mass Flow", colorbar_title="mdot (kg/s)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p2 = heatmap(data.omega_grid, pr_grid, permutedims(data.eta); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Efficiency", colorbar_title="eta (-)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p3 = heatmap(data.omega_grid, pr_grid, permutedims(data.power ./ 1_000.0); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Power", colorbar_title="kW", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p4 = heatmap(data.omega_grid, pr_grid, permutedims(branch_count); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Branch Count", colorbar_title="# branches", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    for plt in (p1, p2, p3, p4)
        _overlay_sampling_points!(plt, overlay_points)
    end
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1300, 950))
    savefig(fig, output_path)
end

function _plot_all(data; output_path::AbstractString, overlay_points=NamedTuple[])
    pt_in = data.diagnostics[1, 1].pt_in
    pr_grid = data.pt_out_grid ./ pt_in
    branch_count = _branch_count_matrix(data.point_results)
    converged = float.(branch_count .> 0.0)
    mdot_min = _candidate_extreme_field_matrix(data.point_results, :mdot, :low)
    mdot_max = _candidate_extreme_field_matrix(data.point_results, :mdot, :high)
    eta_min_branch = _candidate_extreme_field_matrix(data.point_results, :efficiency, :low)
    eta_max_branch = _candidate_extreme_field_matrix(data.point_results, :efficiency, :high)

    p1 = heatmap(data.omega_grid, pr_grid, permutedims(branch_count); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Branch Count", colorbar_title="# branches", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p2 = heatmap(data.omega_grid, pr_grid, permutedims(converged); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Converged Region", colorbar_title="converged", clims=(0.0, 1.0), left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p3 = heatmap(data.omega_grid, pr_grid, permutedims(mdot_min); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Lowest-Flow Branch", colorbar_title="mdot_low (kg/s)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p4 = heatmap(data.omega_grid, pr_grid, permutedims(mdot_max); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Highest-Flow Branch", colorbar_title="mdot_high (kg/s)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p5 = heatmap(data.omega_grid, pr_grid, permutedims(eta_min_branch); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Lowest-Branch Efficiency", colorbar_title="eta_low (-)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p6 = heatmap(data.omega_grid, pr_grid, permutedims(eta_max_branch); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Highest-Branch Efficiency", colorbar_title="eta_high (-)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    for plt in (p1, p2, p3, p4, p5, p6)
        _overlay_sampling_points!(plt, overlay_points)
    end

    fig = plot(p1, p2, p3, p4, p5, p6; layout=(3, 2), size=(1300, 1350))
    savefig(fig, output_path)
end

function main()
    settings = ArgParseSettings(prog="plot_operating_sweep.jl", description="Plot a generic turbomachine operating sweep from a common tabulated performance map.")
    @add_arg_table! settings begin
        "map_path"
            help = "Path to the input common performance-map TOML."
        "--map-group"
            help = "TOML group containing the common performance map."
            arg_type = String
            default = "performance_map"
        "--diagnostics-path"
            help = "Optional path to a performance-map diagnostics TOML on the same grid as the map."
            arg_type = String
        "--diagnostics-group"
            help = "TOML group containing the performance-map diagnostics."
            arg_type = String
            default = "performance_map_diagnostics"
        "--machine-kind"
            help = "Controls default PR bounds and plot labeling: compressor or turbine."
            default = "compressor"
        "--pt-in"
            help = "Inlet stagnation pressure."
            arg_type = Float64
            default = 101325.0
        "--Tt-in"
            help = "Inlet stagnation temperature."
            arg_type = Float64
            default = 288.15
        "--omega-min"
            help = "Minimum shaft speed."
            arg_type = Float64
        "--omega-max"
            help = "Maximum shaft speed."
            arg_type = Float64
        "--n-omega"
            help = "Number of omega samples."
            arg_type = Int
            default = 81
        "--pr-min"
            help = "Minimum pressure ratio pt_out/pt_in."
            arg_type = Float64
        "--pr-max"
            help = "Maximum pressure ratio pt_out/pt_in."
            arg_type = Float64
        "--n-pr"
            help = "Number of PR samples."
            arg_type = Int
            default = 81
        "--branch-selection"
            help = "Branch selection policy for single-branch mode: low, high, nearest."
            default = "nearest"
        "--initial-branch-coordinate"
            help = "Initial branch coordinate used by policy=:nearest."
            arg_type = Float64
        "--track-all-branches"
            help = "Retain and track all branches instead of selecting one."
            action = :store_true
        "--max-match-cost"
            help = "Maximum branch matching cost for all-branch mode."
            arg_type = Float64
            default = Inf
        "--csv"
            help = "Optional CSV output path."
            arg_type = String
        "--output"
            help = "Output plot path."
            arg_type = String
            default = "operating_sweep.png"
    end

    parsed = parse_args(settings)
    machine_kind = Symbol(_parsed_opt(parsed, "machine_kind", "machine-kind"))
    machine_kind in (:compressor, :turbine) || error("machine-kind must be compressor or turbine")

    map = _load_map(parsed["map_path"]; group=_parsed_opt(parsed, "map_group", "map-group", default="performance_map"))
    overlay_points = let path = _parsed_opt(parsed, "diagnostics_path", "diagnostics-path", default=nothing)
        if isnothing(path)
            NamedTuple[]
        else
            diagnostics = _load_diagnostics(String(path); group=_parsed_opt(parsed, "diagnostics_group", "diagnostics-group", default="performance_map_diagnostics"))
            _diagnostic_overlay_points(diagnostics)
        end
    end
    eos = Fluids.ideal_EOS()[:air]
    pt_in = Float64(_parsed_opt(parsed, "pt_in", "pt-in"))
    Tt_in = Float64(_parsed_opt(parsed, "Tt_in", "Tt-in"))
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)

    omega_bounds = _default_omega_bounds(map, Tt_in, pt_in)
    pr_bounds = _default_pr_bounds(map, machine_kind)
    omega_grid = _resolve_grid(_parsed_opt(parsed, "omega_min", "omega-min", default=nothing), _parsed_opt(parsed, "omega_max", "omega-max", default=nothing), Int(_parsed_opt(parsed, "n_omega", "n-omega")), omega_bounds)
    pr_grid = _resolve_grid(_parsed_opt(parsed, "pr_min", "pr-min", default=nothing), _parsed_opt(parsed, "pr_max", "pr-max", default=nothing), Int(_parsed_opt(parsed, "n_pr", "n-pr")), pr_bounds)
    pt_out_grid = pt_in .* pr_grid

    data = TM.solve_operating_sweep(
        map,
        eos;
        pt_in=pt_in,
        ht_in=ht_in,
        omega_grid=omega_grid,
        pt_out_grid=pt_out_grid,
        branch_selection=Symbol(_parsed_opt(parsed, "branch_selection", "branch-selection")),
        initial_branch_coordinate=_parsed_opt(parsed, "initial_branch_coordinate", "initial-branch-coordinate", default=nothing),
        track_all_branches=Bool(_parsed_opt(parsed, "track_all_branches", "track-all-branches", default=false)),
        max_match_cost=Float64(_parsed_opt(parsed, "max_match_cost", "max-match-cost")),
    )

    rows = data.mode == :single ? _single_rows(data) : _all_rows(data)
    csv_path = _parsed_opt(parsed, "csv", "csv", default=nothing)
    if !isnothing(csv_path)
        _write_csv(csv_path, rows)
        println("Wrote CSV to $(csv_path)")
    end

    output_path = _parsed_opt(parsed, "output", "output")
    if data.mode == :single
        _plot_single(data; output_path=output_path, overlay_points=overlay_points)
    else
        _plot_all(data; output_path=output_path, overlay_points=overlay_points)
    end
    println("Saved operating-sweep plot to $(output_path)")
end

main()
