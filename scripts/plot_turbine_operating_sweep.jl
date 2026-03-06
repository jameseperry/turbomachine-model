#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots

const Fluids = TurboMachineModel.Physics.Fluids
const TT = TurboMachineModel.Physics.Turbomachine.Turbine
const U = TurboMachineModel.Utility

function _infer_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".toml"
        return :toml
    end
    error("unsupported map extension $(ext) for path $(path); expected .toml")
end

function _load_turbine_map(path::AbstractString; group::AbstractString="turbine_map")
    _infer_format(path)
    return TT.read_toml(TT.TabulatedTurbinePerformanceMap, path; group=group)
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

_omega_corr_from_omega(map::TT.TabulatedTurbinePerformanceMap, omega::Float64, Tt_in::Float64) =
    TT.corrected_speed(omega, Tt_in, map)

_omega_from_omega_corr(map::TT.TabulatedTurbinePerformanceMap, omega_corr::Float64, Tt_in::Float64) =
    omega_corr * sqrt(Tt_in / map.Tt_ref)

function _state_from_pt_out(
    map::TT.TabulatedTurbinePerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    pt_out::Float64,
    Tt_in::Float64,
)
    vals = TT.turbine_performance_map_from_stagnation(map, omega, pt_in, pt_out, Tt_in)
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    ht_out = ht_in - vals.eta * (ht_in - h2s)
    shaft_power = vals.mdot * (ht_out - ht_in)
    return (
        pt_out=pt_out,
        PR_turb=vals.PR_turb,
        eta=vals.eta,
        mdot=vals.mdot,
        ht_out=ht_out,
        shaft_power=shaft_power,
        valid=(hasproperty(vals, :valid) ? vals.valid : true),
    )
end

function _turbine_pr_bounds(
    map::TT.TabulatedTurbinePerformanceMap,
    omega::Float64,
    Tt_in::Float64,
)
    omega_corr = _omega_corr_from_omega(map, omega, Tt_in)
    domain = TT.performance_map_domain(map)
    if hasproperty(domain, :pr_turb_range)
        pr_lo = domain.pr_turb_range.min(omega_corr)
        pr_hi = domain.pr_turb_range.max(omega_corr)
        return (min(pr_lo, pr_hi), max(pr_lo, pr_hi))
    end
    pr_lo, pr_hi = domain.pr_turb
    return (min(pr_lo, pr_hi), max(pr_lo, pr_hi))
end

function _turbine_torque_roots(
    map::TT.TabulatedTurbinePerformanceMap,
    eos::Fluids.AbstractEOS;
    omega::Float64,
    pt_in::Float64,
    ht_in::Float64,
    Tt_in::Float64,
    tau_load::Float64,
    n_scan::Int=401,
    root_tol::Float64=1e-8,
    prior_roots::AbstractVector{<:Real}=Float64[],
)
    pr_min, pr_max = _turbine_pr_bounds(map, omega, Tt_in)
    pr_min > 0 || return NamedTuple[]
    pr_max > pr_min || return NamedTuple[]
    pt_out_lo = pt_in / pr_max
    pt_out_hi = pt_in / pr_min
    pt_out_hi > pt_out_lo || return NamedTuple[]

    mismatch = function (pt_out)
        state = _state_from_pt_out(map, eos, pt_in, ht_in, omega, pt_out, Tt_in)
        state.valid || return NaN
        return tau_load * omega - state.shaft_power
    end

    roots_pt = U.bracket_bisect_roots(
        mismatch,
        (pt_out_lo, pt_out_hi);
        n_scan=n_scan,
        root_tol=root_tol,
        prior_roots=prior_roots,
    )

    roots = NamedTuple[]
    for pt_out in roots_pt
        state = _state_from_pt_out(map, eos, pt_in, ht_in, omega, pt_out, Tt_in)
        state.valid || continue
        push!(
            roots,
            (
                pt_out=pt_out,
                PR_turb=state.PR_turb,
                eta=state.eta,
                mdot=state.mdot,
                shaft_power=state.shaft_power,
                power_out=-state.shaft_power,
            ),
        )
    end
    sort!(roots; by=r -> r.mdot)
    return roots
end

function _default_tau_load(
    map::TT.TabulatedTurbinePerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Float64,
    Tt_in::Float64,
)
    domain = TT.performance_map_domain(map)
    ωcorr_lo, ωcorr_hi = domain.omega_corr
    ωcorr_hi > ωcorr_lo || return -500.0
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)
    torque_samples = Float64[]

    for ωcorr in range(ωcorr_lo, ωcorr_hi, length=7)
        ω = _omega_from_omega_corr(map, ωcorr, Tt_in)
        ω > 0 || continue
        pr_lo, pr_hi = _turbine_pr_bounds(map, ω, Tt_in)
        pr_hi > pr_lo > 0 || continue
        for pr in range(pr_lo, pr_hi, length=7)
            pt_out = pt_in / pr
            state = _state_from_pt_out(map, eos, pt_in, ht_in, ω, pt_out, Tt_in)
            state.valid || continue
            isfinite(state.shaft_power) || continue
            τ = state.shaft_power / ω
            τ < 0 && push!(torque_samples, τ)
        end
    end

    isempty(torque_samples) && return -500.0
    sort!(torque_samples)
    τ_median = torque_samples[cld(length(torque_samples), 2)]
    # Back off from median extraction torque to improve likelihood of feasible roots.
    return 0.5 * τ_median
end

function _select_roots(roots::Vector{NamedTuple}, branch::Symbol)
    isempty(roots) && return NamedTuple[]
    if branch == :low
        return [first(roots)]
    elseif branch == :high
        return [last(roots)]
    else
        return roots
    end
end

function solve_turbine_operating_sweep(
    map::TT.TabulatedTurbinePerformanceMap,
    eos::Fluids.AbstractEOS;
    omega_min::Float64,
    omega_max::Float64,
    n_points::Int,
    tau_load::Float64,
    pt_in::Float64,
    Tt_in::Float64,
    branch::Symbol=:high,
    root_scan::Int=401,
)
    branch in (:low, :high, :all) || error("branch must be one of: low|high|all")
    n_points >= 2 || error("n_points must be >= 2")

    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)
    omegas = collect(range(omega_min, omega_max, length=n_points))

    if branch == :all
        rows = NamedTuple[]
        diagnostics = NamedTuple[]
        prior_roots = Float64[]
        for omega in omegas
            pr_min, pr_max = _turbine_pr_bounds(map, omega, Tt_in)
            roots = _turbine_torque_roots(
                map,
                eos;
                omega=omega,
                pt_in=pt_in,
                ht_in=ht_in,
                Tt_in=Tt_in,
                tau_load=tau_load,
                n_scan=root_scan,
                prior_roots=prior_roots,
            )
            push!(
                diagnostics,
                (
                    omega=omega,
                    converged=!isempty(roots),
                    n_roots=length(roots),
                    pr_turb_min=pr_min,
                    pr_turb_max=pr_max,
                    failure_reason=isempty(roots) ? "no_torque_root_on_pr_domain" : "none",
                ),
            )
            if isempty(roots)
                push!(
                    rows,
                    (
                        omega=omega,
                        branch_id=0,
                        PR_turb=NaN,
                        eta=NaN,
                        mdot=NaN,
                        power_out=NaN,
                        converged=false,
                        pt_out=NaN,
                    ),
                )
                continue
            end
            for (k, r) in enumerate(roots)
                push!(
                    rows,
                    (
                        omega=omega,
                        branch_id=k,
                        PR_turb=r.PR_turb,
                        eta=r.eta,
                        mdot=r.mdot,
                        power_out=r.power_out,
                        converged=true,
                        pt_out=r.pt_out,
                    ),
                )
            end
            prior_roots = [r.pt_out for r in roots]
        end
        return (mode=:all, branch=:all, rows=rows, diagnostics=diagnostics)
    end

    prs = fill(NaN, n_points)
    etas = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    pt_outs = fill(NaN, n_points)
    converged = fill(false, n_points)
    diagnostics = NamedTuple[]
    prior_roots = Float64[]

    for (i, omega) in enumerate(omegas)
        pr_min, pr_max = _turbine_pr_bounds(map, omega, Tt_in)
        roots = _turbine_torque_roots(
            map,
            eos;
            omega=omega,
            pt_in=pt_in,
            ht_in=ht_in,
            Tt_in=Tt_in,
            tau_load=tau_load,
            n_scan=root_scan,
            prior_roots=prior_roots,
        )
        push!(
            diagnostics,
            (
                omega=omega,
                converged=!isempty(roots),
                n_roots=length(roots),
                pr_turb_min=pr_min,
                pr_turb_max=pr_max,
                failure_reason=isempty(roots) ? "no_torque_root_on_pr_domain" : "none",
            ),
        )
        isempty(roots) && continue
        root = only(_select_roots(roots, branch))
        prs[i] = root.PR_turb
        etas[i] = root.eta
        mdots[i] = root.mdot
        powers[i] = root.power_out
        pt_outs[i] = root.pt_out
        converged[i] = true
        prior_roots = [r.pt_out for r in roots]
    end

    return (
        mode=:single,
        branch=branch,
        omegas=omegas,
        prs=prs,
        etas=etas,
        mdots=mdots,
        powers=powers,
        pt_outs=pt_outs,
        converged=converged,
        diagnostics=diagnostics,
    )
end

_diagnostics(data) = hasproperty(data, :diagnostics) ? data.diagnostics : NamedTuple[]

function _plot_turbine_operating_sweep_data(
    data;
    output_path::AbstractString="turbine_operating_sweep.png",
    branch_match_cost::Float64=0.5,
)
    branch_match_cost >= 0 || error("branch_match_cost must be >= 0")

    if data.mode == :single
        omega_axis = data.omegas
        powers_kw = data.powers ./ 1_000.0

        p1 = plot(
            omega_axis,
            data.prs;
            xlabel="shaft speed omega",
            ylabel="expansion ratio Pt_in/Pt_out",
            title="Expansion Ratio vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )
        p2 = plot(
            omega_axis,
            data.mdots;
            xlabel="shaft speed omega",
            ylabel="mass flow rate mdot (kg/s)",
            title="Mass Flow vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )
        p3 = plot(
            omega_axis,
            powers_kw;
            xlabel="shaft speed omega",
            ylabel="power output (kW)",
            title="Power Output vs Shaft Speed",
            lw=2,
            marker=:circle,
            label=false,
        )
        p4 = plot(
            omega_axis,
            data.etas;
            xlabel="shaft speed omega",
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
        println("Converged points: $(count(data.converged)) / $(length(data.omegas))")
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
                PR_turb=Float64(r.PR_turb),
                eta=Float64(r.eta),
                mdot=Float64(r.mdot),
                power_kw=Float64(r.power_out) / 1_000.0,
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
        pr_vals = [r.PR_turb for r in flat_roots]
        eta_vals = [r.eta for r in flat_roots]
        mdot_scale = max(maximum(mdot_vals) - minimum(mdot_vals), 1e-6)
        pr_scale = max(maximum(pr_vals) - minimum(pr_vals), 1e-6)
        eta_scale = max(maximum(eta_vals) - minimum(eta_vals), 1e-6)
        U.track_branches(
            omegas,
            roots_by_condition;
            distance=(a, b) ->
                abs(a.mdot - b.mdot) / mdot_scale +
                abs(a.PR_turb - b.PR_turb) / pr_scale +
                abs(a.eta - b.eta) / eta_scale,
            max_match_cost=branch_match_cost,
        )
    end

    p1 = plot(
        xlabel="shaft speed omega",
        ylabel="expansion ratio Pt_in/Pt_out",
        title="Expansion Ratio vs Shaft Speed",
        label=false,
    )
    p2 = plot(
        xlabel="shaft speed omega",
        ylabel="mass flow rate mdot (kg/s)",
        title="Mass Flow vs Shaft Speed",
        label=false,
    )
    p3 = plot(
        xlabel="shaft speed omega",
        ylabel="power output (kW)",
        title="Power Output vs Shaft Speed",
        label=false,
    )
    p4 = plot(
        xlabel="shaft speed omega",
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
            pr_line[i] = point.root.PR_turb
            mdot_line[i] = point.root.mdot
            power_line[i] = point.root.power_kw
            eta_line[i] = point.root.eta
        end
        plot!(p1, omegas, pr_line; lw=2, marker=:circle, ms=3, label=false)
        plot!(p2, omegas, mdot_line; lw=2, marker=:circle, ms=3, label=false)
        plot!(p3, omegas, power_line; lw=2, marker=:circle, ms=3, label=false)
        plot!(p4, omegas, eta_line; lw=2, marker=:circle, ms=3, label=false)
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
    println("Converged rows: $(count(r -> r.converged, data.rows)) / $(length(data.rows))")
    println("Tracked branches: $(length(branch_tracks.branch_ids))")
    println("Saved operating-point sweep plot to: $output_path")
end

function write_turbine_operating_sweep_csv(
    data;
    output_path::AbstractString="turbine_operating_sweep.csv",
)
    diag_by_omega = Dict{Float64,NamedTuple}()
    for d in _diagnostics(data)
        diag_by_omega[Float64(d.omega)] = d
    end
    open(output_path, "w") do io
        println(
            io,
            "omega,branch_id,PR_turb,eta,mdot,power_out_kw,pt_out_pa,converged,failure_reason,n_roots,pr_turb_min,pr_turb_max",
        )
        if data.mode == :single
            for i in eachindex(data.omegas)
                branch_id = data.branch == :low ? 1 : 2
                d = get(
                    diag_by_omega,
                    Float64(data.omegas[i]),
                    (
                        failure_reason="none",
                        n_roots=0,
                        pr_turb_min=NaN,
                        pr_turb_max=NaN,
                    ),
                )
                println(
                    io,
                    "$(data.omegas[i]),$branch_id,$(data.prs[i]),$(data.etas[i]),$(data.mdots[i]),$(data.powers[i] / 1_000.0),$(data.pt_outs[i]),$(data.converged[i]),$(d.failure_reason),$(d.n_roots),$(d.pr_turb_min),$(d.pr_turb_max)",
                )
            end
        else
            for r in data.rows
                d = get(
                    diag_by_omega,
                    Float64(r.omega),
                    (
                        failure_reason="none",
                        n_roots=0,
                        pr_turb_min=NaN,
                        pr_turb_max=NaN,
                    ),
                )
                println(
                    io,
                    "$(r.omega),$(r.branch_id),$(r.PR_turb),$(r.eta),$(r.mdot),$(r.power_out / 1_000.0),$(r.pt_out),$(r.converged),$(d.failure_reason),$(d.n_roots),$(d.pr_turb_min),$(d.pr_turb_max)",
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
end

function print_turbine_operating_sweep_failures(data)
    failures = filter(d -> !d.converged, _diagnostics(data))
    isempty(failures) && return
    println("Failed points diagnostics:")
    for d in failures
        println(
            "  omega=$(d.omega), reason=$(d.failure_reason), roots=$(d.n_roots), PR_turb_range=[$(d.pr_turb_min), $(d.pr_turb_max)]",
        )
    end
end

function _main(args::Vector{String}=ARGS)
    settings = ArgParseSettings(
        prog="plot_turbine_operating_sweep.jl",
        description="Solve and plot turbine operating sweep from a tabulated performance map file.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input turbine performance map (.toml)"
            required = true
        "--output"
            help = "output plot path"
            arg_type = String
            default = "turbine_operating_sweep.png"
        "--csv"
            help = "optional CSV output path"
            arg_type = String
        "--map-group"
            help = "input map group/table"
            arg_type = String
            default = "turbine_map"
        "--tau-load"
            help = "load torque (N*m), typically negative for power extraction; default is estimated from map/inlet state"
            arg_type = Float64
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
        "--branch"
            help = "root branch to track: low, high, or all"
            arg_type = String
            default = "high"
        "--root-scan"
            help = "number of scan points used to bracket pt_out roots per speed"
            arg_type = Int
            default = 401
        "--branch-match-cost"
            help = "max normalized match cost when tracking branches for --branch all plots"
            arg_type = Float64
            default = 0.5
    end

    parsed = parse_args(args, settings)
    map_group = something(_parsed_opt(parsed, "map_group", "map-group"), "turbine_map")
    csv_output = _parsed_opt(parsed, "csv", "csv")
    branch = _parse_branch(something(_parsed_opt(parsed, "branch", "branch"), "high"))
    output = something(_parsed_opt(parsed, "output", "output"), "turbine_operating_sweep.png")
    tau_load_arg = _parsed_opt(parsed, "tau_load", "tau-load")
    omega_min_arg = _parsed_opt(parsed, "omega_min", "omega-min")
    omega_max_arg = _parsed_opt(parsed, "omega_max", "omega-max")
    n_points = Int(something(_parsed_opt(parsed, "n_points", "n-points"), 25))
    pt_in = Float64(something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0))
    tt_in = Float64(something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15))
    root_scan = Int(something(_parsed_opt(parsed, "root_scan", "root-scan"), 401))
    branch_match_cost = Float64(something(_parsed_opt(parsed, "branch_match_cost", "branch-match-cost"), 0.5))

    map = _load_turbine_map(parsed["map_path"]; group=map_group)
    domain = TT.performance_map_domain(map)
    ωcorr_min, ωcorr_max = domain.omega_corr
    omega_default_min = _omega_from_omega_corr(map, ωcorr_min, tt_in)
    omega_default_max = _omega_from_omega_corr(map, ωcorr_max, tt_in)
    omega_min = Float64(something(omega_min_arg, omega_default_min))
    omega_max = Float64(something(omega_max_arg, omega_default_max))
    tau_load = if isnothing(tau_load_arg)
        _default_tau_load(map, Fluids.ideal_EOS()[:air]; pt_in=pt_in, Tt_in=tt_in)
    else
        Float64(tau_load_arg)
    end
    println("Using omega sweep range [$(omega_min), $(omega_max)] rad/s and tau_load=$(tau_load) N*m")

    data = solve_turbine_operating_sweep(
        map,
        Fluids.ideal_EOS()[:air];
        omega_min=omega_min,
        omega_max=omega_max,
        n_points=n_points,
        tau_load=tau_load,
        pt_in=pt_in,
        Tt_in=tt_in,
        branch=branch,
        root_scan=root_scan,
    )

    _plot_turbine_operating_sweep_data(
        data;
        output_path=output,
        branch_match_cost=branch_match_cost,
    )

    if !isnothing(csv_output)
        write_turbine_operating_sweep_csv(data; output_path=String(csv_output))
    end
    print_turbine_operating_sweep_failures(data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
