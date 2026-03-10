#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots
using Statistics: mean
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

function _plot_single(data; output_path::AbstractString)
    pt_in = data.diagnostics[1, 1].pt_in
    pr_grid = data.pt_out_grid ./ pt_in
    branch_count = _branch_count_matrix(data.point_results)

    p1 = heatmap(data.omega_grid, pr_grid, permutedims(data.mdot); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Mass Flow", colorbar_title="mdot (kg/s)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p2 = heatmap(data.omega_grid, pr_grid, permutedims(data.eta); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Efficiency", colorbar_title="eta (-)", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p3 = heatmap(data.omega_grid, pr_grid, permutedims(data.power ./ 1_000.0); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Power", colorbar_title="kW", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    p4 = heatmap(data.omega_grid, pr_grid, permutedims(branch_count); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Branch Count", colorbar_title="# branches", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1300, 950))
    savefig(fig, output_path)
end

function _plot_all(data; output_path::AbstractString)
    pt_in = data.diagnostics[1, 1].pt_in
    pr_grid = data.pt_out_grid ./ pt_in
    branch_count = _branch_count_matrix(data.point_results)
    good_rows = filter(r -> r.converged, data.rows)

    p1 = heatmap(data.omega_grid, pr_grid, permutedims(branch_count); xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Branch Count", colorbar_title="# branches", left_margin=8mm, right_margin=10mm, bottom_margin=8mm, top_margin=8mm)

    if !isempty(good_rows)
        p2 = _plot_branch_field(good_rows, :mdot; xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Mass Flow Branches", colorbar_title="mdot (kg/s)")
        p3 = _plot_branch_field(good_rows, :eta; xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Efficiency Branches", colorbar_title="eta (-)")
        p4 = _plot_branch_field(good_rows, :power_kw; xlabel="omega (rad/s)", ylabel="PR = pt_out / pt_in", title="Power Branches", colorbar_title="kW")
    else
        p2 = plot(title="Mass Flow Branches", axis=false)
        p3 = plot(title="Efficiency Branches", axis=false)
        p4 = plot(title="Power Branches", axis=false)
    end

    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1300, 950))
    savefig(fig, output_path)
end

function _plot_branch_field(
    rows::Vector{<:NamedTuple},
    field::Symbol;
    xlabel::AbstractString,
    ylabel::AbstractString,
    title::AbstractString,
    colorbar_title::AbstractString,
)
    branch_ids = sort!(unique([r.branch_id for r in rows if r.branch_id > 0]))
    p = plot(
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        colorbar_title=colorbar_title,
        legend=false,
        left_margin=8mm,
        right_margin=10mm,
        bottom_margin=8mm,
        top_margin=8mm,
    )

    all_vals = Float64[]
    for row in rows
        value = field === :power_kw ? row.power / 1_000.0 : getproperty(row, field)
        isfinite(value) && push!(all_vals, Float64(value))
    end
    isempty(all_vals) && return p

    vmin = minimum(all_vals)
    vmax = maximum(all_vals)
    span = vmax - vmin
    color_for(value) = cgrad(:viridis)[span <= 0 ? 0.5 : clamp((value - vmin) / span, 0.0, 1.0)]

    first_series = true
    for branch_id in branch_ids
        branch_rows = sort(filter(r -> r.branch_id == branch_id, rows), by=r -> r.omega)
        length(branch_rows) >= 2 || continue

        omega = [r.omega for r in branch_rows]
        pr = [r.PR for r in branch_rows]
        values = [field === :power_kw ? r.power / 1_000.0 : getproperty(r, field) for r in branch_rows]
        mean_value = mean(values)

        plot!(
            p,
            omega,
            pr;
            color=color_for(mean_value),
            lw=2,
            label=false,
        )

        if first_series
            scatter!(
                p,
                [omega[1]],
                [pr[1]];
                zcolor=[mean_value],
                color=:viridis,
                ms=0,
                markerstrokewidth=0,
                label=false,
                colorbar=true,
                colorbar_title=colorbar_title,
                clims=(vmin, vmax),
            )
            first_series = false
        end
    end

    return p
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
        _plot_single(data; output_path=output_path)
    else
        _plot_all(data; output_path=output_path)
    end
    println("Saved operating-sweep plot to $(output_path)")
end

main()
