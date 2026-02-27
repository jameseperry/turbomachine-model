#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots

const TM = TurboMachineModel.Physics.Turbomachine.Compressor

function _infer_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".toml"
        return :toml
    end
    error("unsupported map extension $(ext) for path $(path); expected .toml")
end

function _load_compressor_map(
    path::AbstractString;
    group::Union{Nothing,AbstractString}=nothing,
)
    _infer_format(path)
    groups = isnothing(group) ? ("compressor_map", "compressor_analytic_map") : (group,)
    for g in groups
        try
            return TM.read_toml(TM.TabulatedCompressorPerformanceMap, path; group=g)
        catch
        end
        try
            return TM.read_toml(TM.AnalyticCompressorPerformanceMap, path; group=g)
        catch
        end
    end
    error(
        "failed to load compressor map from $(path); expected tabulated group compressor_map or analytic group compressor_analytic_map",
    )
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _parse_branch(raw::AbstractString)
    b = Symbol(lowercase(strip(raw)))
    b in (:low, :high, :all) || error("branch must be one of: low|high|all")
    return b
end

function _plot_compressor_operating_sweep_data(
    data;
    output_path::AbstractString="compressor_operating_sweep.png",
    branch_match_cost::Float64=0.5,
)
    branch_match_cost >= 0 || error("branch_match_cost must be >= 0")

    if data.mode == :single
        omega_rad_s = data.omegas

        p1 = plot(
            omega_rad_s,
            data.prs;
            xlabel="shaft speed omega (rad/s)",
            ylabel="compression ratio Pt_out/Pt_in",
            title="Compression Ratio vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )

        p2 = plot(
            omega_rad_s,
            data.mdots;
            xlabel="shaft speed omega (rad/s)",
            ylabel="mass flow rate mdot (kg/s)",
            title="Mass Flow vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )

        powers_kw = data.powers ./ 1_000.0

        p3 = plot(
            omega_rad_s,
            powers_kw;
            xlabel="shaft speed omega (rad/s)",
            ylabel="power consumption (kW)",
            title="Power Consumption vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )

        p4 = plot(
            omega_rad_s,
            data.etas;
            xlabel="shaft speed omega (rad/s)",
            ylabel="adiabatic efficiency eta (-)",
            title="Adiabatic Efficiency vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )

        fig = plot(
            p1,
            p2,
            p3,
            p4;
            layout=(2, 2),
            size=(1200, 900),
            left_margin=8Plots.mm,
            right_margin=8Plots.mm,
            top_margin=8Plots.mm,
            bottom_margin=8Plots.mm,
        )
        savefig(fig, output_path)

        n_converged = count(data.converged)
        println("Converged points: $n_converged / $(length(data.omegas))")
        println("Saved operating-point sweep plot to: $output_path")
        return
    end

    data.mode == :all || error("unsupported sweep mode: $(data.mode)")

    omegas = sort!(collect(Set(Float64(r.omega) for r in data.rows)))
    n_omega = length(omegas)
    omega_to_idx = Dict(omega => i for (i, omega) in enumerate(omegas))
    roots_by_condition = [NamedTuple[] for _ in 1:n_omega]
    for r in data.rows
        r.converged || continue
        i = omega_to_idx[Float64(r.omega)]
        push!(
            roots_by_condition[i],
            (
                PR=Float64(r.PR),
                eta=Float64(r.eta),
                mdot=Float64(r.mdot),
                power_kw=Float64(r.power) / 1_000.0,
            ),
        )
    end

    flat_roots = NamedTuple[]
    for roots in roots_by_condition
        append!(flat_roots, roots)
    end

    branch_tracks = if isempty(flat_roots)
        (branch_ids=Int[], branches=Dict{Int,Vector{NamedTuple}}())
    else
        mdot_vals = [r.mdot for r in flat_roots]
        pr_vals = [r.PR for r in flat_roots]
        eta_vals = [r.eta for r in flat_roots]
        mdot_scale = max(maximum(mdot_vals) - minimum(mdot_vals), 1e-6)
        pr_scale = max(maximum(pr_vals) - minimum(pr_vals), 1e-6)
        eta_scale = max(maximum(eta_vals) - minimum(eta_vals), 1e-6)

        TurboMachineModel.Utility.track_branches(
            omegas,
            roots_by_condition;
            distance=(a, b) ->
                abs(a.mdot - b.mdot) / mdot_scale +
                abs(a.PR - b.PR) / pr_scale +
                abs(a.eta - b.eta) / eta_scale,
            max_match_cost=branch_match_cost,
        )
    end

    omega_rad_s = omegas
    p1 = plot(
        xlabel="shaft speed omega (rad/s)",
        ylabel="compression ratio Pt_out/Pt_in",
        title="Compression Ratio vs Shaft Speed",
        label=false,
    )
    p2 = plot(
        xlabel="shaft speed omega (rad/s)",
        ylabel="mass flow rate mdot (kg/s)",
        title="Mass Flow vs Shaft Speed",
        label=false,
    )
    p3 = plot(
        xlabel="shaft speed omega (rad/s)",
        ylabel="power consumption (kW)",
        title="Power Consumption vs Shaft Speed",
        label=false,
    )
    p4 = plot(
        xlabel="shaft speed omega (rad/s)",
        ylabel="adiabatic efficiency eta (-)",
        title="Adiabatic Efficiency vs Shaft Speed",
        label=false,
    )

    for bid in branch_tracks.branch_ids
        pr_line = fill(NaN, n_omega)
        mdot_line = fill(NaN, n_omega)
        power_line = fill(NaN, n_omega)
        eta_line = fill(NaN, n_omega)
        for point in branch_tracks.branches[bid]
            i = point.condition_idx
            pr_line[i] = point.root.PR
            mdot_line[i] = point.root.mdot
            power_line[i] = point.root.power_kw
            eta_line[i] = point.root.eta
        end
        plot!(p1, omega_rad_s, pr_line; lw=2, marker=:circle, ms=3, label=false)
        plot!(p2, omega_rad_s, mdot_line; lw=2, marker=:circle, ms=3, label=false)
        plot!(p3, omega_rad_s, power_line; lw=2, marker=:circle, ms=3, label=false)
        plot!(p4, omega_rad_s, eta_line; lw=2, marker=:circle, ms=3, label=false)
    end

    fig = plot(
        p1,
        p2,
        p3,
        p4;
        layout=(2, 2),
        size=(1200, 900),
        left_margin=8Plots.mm,
        right_margin=8Plots.mm,
        top_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    savefig(fig, output_path)

    n_converged = count(r -> r.converged, data.rows)
    println("Converged rows: $n_converged / $(length(data.rows))")
    println("Tracked branches: $(length(branch_tracks.branch_ids))")
    println("Branch match cost: $branch_match_cost")
    println("Saved operating-point sweep plot to: $output_path")
end

function write_compressor_operating_sweep_csv(
    data;
    output_path::AbstractString="compressor_operating_sweep.csv",
)
    open(output_path, "w") do io
        println(io, "omega_rad_per_s,branch_id,PR,eta,mdot,power_kw,converged,backoff_used")
        if data.mode == :single
            for i in eachindex(data.omegas)
                branch_id = data.branch == :low ? 1 : 2
                power_kw = data.powers[i] / 1_000.0
                println(
                    io,
                    "$(data.omegas[i]),$branch_id,$(data.prs[i]),$(data.etas[i]),$(data.mdots[i]),$(power_kw),$(data.converged[i]),$(data.backoff_used[i])",
                )
            end
        else
            for r in data.rows
                power_kw = r.power / 1_000.0
                println(
                    io,
                    "$(r.omega),$(r.branch_id),$(r.PR),$(r.eta),$(r.mdot),$(power_kw),$(r.converged),$(r.backoff_used)",
                )
            end
        end
    end

    n_converged = if data.mode == :single
        count(data.converged)
    else
        count(r -> r.converged, data.rows)
    end
    n_total = if data.mode == :single
        length(data.omegas)
    else
        length(data.rows)
    end

    println("Converged rows: $n_converged / $n_total")
    println("Saved operating-point sweep CSV to: $output_path")
    return output_path
end

function _main(args::Vector{String}=ARGS)
    settings = ArgParseSettings(
        prog="plot_compressor_operating_sweep.jl",
        description="Solve and plot compressor operating sweep from a compressor performance map file.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input compressor performance map (.toml)"
            required = true
        "--output"
            help = "output plot/csv path"
            arg_type = String
        "--csv"
            help = "write CSV output instead of generating a plot"
            action = :store_true
        "--map-group"
            help = "input map group/table"
            arg_type = String
        "--target-pr"
            help = "target compressor pressure ratio Pt_out/Pt_in"
            arg_type = Float64
            default = 2.0
        "--branch"
            help = "root branch to track: low, high, or all"
            arg_type = String
            default = "high"
        "--disable-pr-backoff"
            help = "disable pt_out backoff when target PR has no root"
            action = :store_true
        "--backoff-min-pt-out"
            help = "minimum pt_out to consider in backoff search (Pa); default=pt_in"
            arg_type = Float64
        "--backoff-max-pt-out"
            help = "maximum pt_out to consider in backoff search (Pa); default=pt_in*target_pr"
            arg_type = Float64
        "--backoff-pt-out-tol"
            help = "pt_out tolerance for backoff bisection (Pa)"
            arg_type = Float64
            default = 50.0
        "--backoff-max-iters"
            help = "maximum bisection iterations for backoff"
            arg_type = Int
            default = 24
        "--omega-min"
            help = "minimum shaft speed omega (rad/s); default comes from map domain"
            arg_type = Float64
        "--omega-max"
            help = "maximum shaft speed omega (rad/s); default comes from map domain"
            arg_type = Float64
        "--n-points"
            help = "number of sweep points"
            arg_type = Int
            default = 25
        "--pt-in"
            help = "inlet total pressure (Pa)"
            arg_type = Float64
            default = 101_325.0
        "--tt-in"
            help = "inlet total temperature (K)"
            arg_type = Float64
            default = 288.15
        "--branch-match-cost"
            help = "max normalized match cost when tracking branches for --branch all plots"
            arg_type = Float64
            default = 0.5
    end

    parsed = parse_args(args, settings)

    csv_mode = something(_parsed_opt(parsed, "csv", "csv"), false)
    map_group = _parsed_opt(parsed, "map_group", "map-group")
    branch = _parse_branch(something(_parsed_opt(parsed, "branch", "branch"), "high"))

    output_default = csv_mode ? "compressor_operating_sweep.csv" : "compressor_operating_sweep.png"
    output = something(_parsed_opt(parsed, "output", "output"), output_default)

    target_pr = something(_parsed_opt(parsed, "target_pr", "target-pr"), 2.0)
    pr_backoff = !something(_parsed_opt(parsed, "disable_pr_backoff", "disable-pr-backoff"), false)
    backoff_min_pt_out = _parsed_opt(parsed, "backoff_min_pt_out", "backoff-min-pt-out")
    backoff_max_pt_out = _parsed_opt(parsed, "backoff_max_pt_out", "backoff-max-pt-out")
    backoff_pt_out_tol = something(_parsed_opt(parsed, "backoff_pt_out_tol", "backoff-pt-out-tol"), 50.0)
    backoff_max_iters = something(_parsed_opt(parsed, "backoff_max_iters", "backoff-max-iters"), 24)

    omega_min_arg = _parsed_opt(parsed, "omega_min", "omega-min")
    omega_max_arg = _parsed_opt(parsed, "omega_max", "omega-max")
    n_points = something(_parsed_opt(parsed, "n_points", "n-points"), 25)
    pt_in = something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0)
    tt_in = something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15)
    branch_match_cost = something(_parsed_opt(parsed, "branch_match_cost", "branch-match-cost"), 0.5)

    map = _load_compressor_map(parsed["map_path"]; group=map_group)
    domain = TM.performance_map_domain(map)
    omega_default_min, omega_default_max = domain.omega_corr
    omega_min = something(omega_min_arg, omega_default_min)
    omega_max = something(omega_max_arg, omega_default_max)
    data = TM.solve_compressor_operating_sweep(
        map;
        omega_min=omega_min,
        omega_max=omega_max,
        n_points=n_points,
        pt_in=pt_in,
        Tt_in=tt_in,
        target_pr=target_pr,
        branch=branch,
        pr_backoff=pr_backoff,
        backoff_min_pt_out=backoff_min_pt_out,
        backoff_max_pt_out=backoff_max_pt_out,
        backoff_pt_out_tol=backoff_pt_out_tol,
        backoff_max_iters=backoff_max_iters,
    )

    if csv_mode
        write_compressor_operating_sweep_csv(data; output_path=output)
    else
        _plot_compressor_operating_sweep_data(
            data;
            output_path=output,
            branch_match_cost=branch_match_cost,
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
