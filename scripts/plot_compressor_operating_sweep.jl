#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots

const Fluids = TurboMachineModel.Physics.Fluids
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

function _inner_operating_point(
    map::TM.AbstractCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    pt_out::Float64,
    mdot_guess::Float64,
    ht_out_guess::Float64,
    tau_guess::Float64,
)
    return TM.solve_compressor_operating_point(
        map,
        eos;
        pt_in=pt_in,
        ht_in=ht_in,
        pt_out=pt_out,
        omega=omega,
        mdot_guess=mdot_guess,
        ht_out_guess=ht_out_guess,
        tau_guess=tau_guess,
    )
end

function _solve_at_pt_out(
    map::TM.AbstractCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    Tt_in::Float64,
    omega::Float64,
    pt_out::Float64,
    seed::Union{Nothing,NamedTuple}=nothing,
)
    target_pr = pt_out / pt_in
    branch_seed = _map_seed_for_target_pr(
        map,
        eos,
        pt_in,
        ht_in,
        Tt_in,
        omega,
        target_pr,
        isnothing(seed) ? nothing : seed.mdot,
    )

    if !isnothing(seed)
        sol = _inner_operating_point(
            map,
            eos,
            pt_in,
            ht_in,
            omega,
            pt_out,
            seed.mdot,
            seed.ht_out,
            seed.tau,
        )
        sol.converged && return (sol=sol, branch_seed=branch_seed)
    end

    sol = _inner_operating_point(
        map,
        eos,
        pt_in,
        ht_in,
        omega,
        pt_out,
        branch_seed.mdot,
        branch_seed.ht_out,
        branch_seed.tau,
    )
    return (sol=sol, branch_seed=branch_seed)
end

function _map_seed_for_target_pr(
    map::TM.AbstractCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    Tt_in::Float64,
    omega::Float64,
    target_pr::Float64,
    preferred_mdot::Union{Nothing,Float64}=nothing,
)
    domain = TM.performance_map_domain(map)
    omega_corr = TM.corrected_speed(omega, Tt_in, map)
    flow_range = domain.mdot_corr_flow_range
    m_surge = flow_range.surge(omega_corr)
    m_choke = flow_range.choke(omega_corr)
    m_lo = min(m_surge, m_choke)
    m_hi = max(m_surge, m_choke)

    m_grid = collect(range(m_lo, m_hi, length=41))
    pr_vals = similar(m_grid)
    eta_vals = similar(m_grid)
    for i in eachindex(m_grid)
        vals = TM.compressor_performance_map(map, omega_corr, m_grid[i])
        pr_vals[i] = vals.PR
        eta_vals[i] = vals.eta
    end

    # Find mdot_corr roots where PR(mdot_corr) = target_pr.
    # If multiple roots exist, prefer continuity from prior mdot; otherwise use
    # the higher-flow (choke-side) branch.
    mdot_corr_candidates = Float64[]
    eta_candidates = Float64[]
    for i in 1:(length(m_grid) - 1)
        f1 = pr_vals[i] - target_pr
        f2 = pr_vals[i + 1] - target_pr
        if f1 == 0.0
            push!(mdot_corr_candidates, m_grid[i])
            push!(eta_candidates, eta_vals[i])
            continue
        end
        if f1 * f2 > 0.0
            continue
        end
        a = m_grid[i]
        b = m_grid[i + 1]
        fa = f1
        mdot_root = 0.5 * (a + b)
        eta_root = 0.5 * (eta_vals[i] + eta_vals[i + 1])
        for _ in 1:40
            mid = 0.5 * (a + b)
            vals = TM.compressor_performance_map(map, omega_corr, mid)
            fm = vals.PR - target_pr
            eta_root = vals.eta
            mdot_root = mid
            if abs(fm) <= 1e-8
                break
            end
            if fa * fm <= 0.0
                b = mid
            else
                a = mid
                fa = fm
            end
        end
        push!(mdot_corr_candidates, mdot_root)
        push!(eta_candidates, eta_root)
    end

    mdot_corr_guess = m_grid[argmin(abs.(pr_vals .- target_pr))]
    eta_guess = eta_vals[argmin(abs.(pr_vals .- target_pr))]
    if !isempty(mdot_corr_candidates)
        if preferred_mdot !== nothing
            preferred_mdot_corr = TM.corrected_flow(preferred_mdot, Tt_in, pt_in, map)
            k = argmin(abs.(mdot_corr_candidates .- preferred_mdot_corr))
            mdot_corr_guess = mdot_corr_candidates[k]
            eta_guess = eta_candidates[k]
        else
            k = argmax(mdot_corr_candidates)
            mdot_corr_guess = mdot_corr_candidates[k]
            eta_guess = eta_candidates[k]
        end
    end

    corr_to_phys_scale = (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
    mdot_guess = mdot_corr_guess * corr_to_phys_scale
    eta_safe = max(eta_guess, 1e-3)
    pt_out = pt_in * target_pr
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    ht_out_guess = ht_in + (h2s - ht_in) / eta_safe
    tau_guess = mdot_guess * (ht_out_guess - ht_in) / max(omega, 1e-6)

    return (mdot=mdot_guess, ht_out=ht_out_guess, tau=tau_guess)
end

function _solve_with_pr_backoff(
    map::TM.AbstractCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    Tt_in::Float64,
    omega::Float64,
    target_pt_out::Float64,
    seed::Union{Nothing,NamedTuple};
    min_pt_out::Float64,
    max_pt_out::Float64,
    pt_out_tol::Float64=50.0,
    max_iters::Int=24,
)
    pt_out_tol > 0 || error("pt_out_tol must be > 0")
    max_iters >= 1 || error("max_iters must be >= 1")

    lo = max(min_pt_out, pt_in)
    hi = min(max_pt_out, target_pt_out)
    lo <= hi || return (converged=false, reason="invalid backoff range")

    hi_try = _solve_at_pt_out(map, eos, pt_in, ht_in, Tt_in, omega, hi, seed)
    if hi_try.sol.converged
        return (converged=true, sol=hi_try.sol, pt_out=hi, used_backoff=(hi < target_pt_out))
    end

    lo_try = _solve_at_pt_out(map, eos, pt_in, ht_in, Tt_in, omega, lo, seed)
    if !lo_try.sol.converged
        return (converged=false, reason="no converged solution in backoff range")
    end

    lo_sol = lo_try.sol
    lo_seed = (mdot=lo_sol.mdot, ht_out=lo_sol.ht_out, tau=lo_sol.tau)

    for _ in 1:max_iters
        if (hi - lo) <= pt_out_tol
            break
        end
        mid = 0.5 * (lo + hi)
        mid_try = _solve_at_pt_out(map, eos, pt_in, ht_in, Tt_in, omega, mid, lo_seed)
        if mid_try.sol.converged
            lo = mid
            lo_sol = mid_try.sol
            lo_seed = (mdot=lo_sol.mdot, ht_out=lo_sol.ht_out, tau=lo_sol.tau)
        else
            hi = mid
        end
    end

    return (converged=true, sol=lo_sol, pt_out=lo, used_backoff=true)
end

function sweep_compressor_operating_points(
    map::TM.AbstractCompressorPerformanceMap;
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
    target_pr::Float64=2.0,
    pr_backoff::Bool=true,
    backoff_min_pt_out::Union{Nothing,Float64}=nothing,
    backoff_max_pt_out::Union{Nothing,Float64}=nothing,
    backoff_pt_out_tol::Float64=50.0,
    backoff_max_iters::Int=24,
)
    target_pr > 1.0 || error("target_pr must be > 1.0")
    eos = Fluids.ideal_EOS()[:air]
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)
    target_pt_out = pt_in * target_pr
    min_pt_out = isnothing(backoff_min_pt_out) ? pt_in : backoff_min_pt_out
    max_pt_out = isnothing(backoff_max_pt_out) ? target_pt_out : backoff_max_pt_out

    omegas = collect(range(omega_min, omega_max, length=n_points))
    prs = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    converged = fill(false, n_points)
    backoff_used = fill(false, n_points)

    seed::Union{Nothing,NamedTuple{(:mdot, :ht_out, :tau),Tuple{Float64,Float64,Float64}}} = nothing

    for (i, omega) in enumerate(omegas)
        direct_try = _solve_at_pt_out(
            map,
            eos,
            pt_in,
            ht_in,
            Tt_in,
            omega,
            target_pt_out,
            seed,
        )
        sol = direct_try.sol

        if !sol.converged && pr_backoff
            backoff = _solve_with_pr_backoff(
                map,
                eos,
                pt_in,
                ht_in,
                Tt_in,
                omega,
                target_pt_out,
                seed;
                min_pt_out=min_pt_out,
                max_pt_out=max_pt_out,
                pt_out_tol=backoff_pt_out_tol,
                max_iters=backoff_max_iters,
            )
            if backoff.converged
                sol = backoff.sol
                backoff_used[i] = backoff.used_backoff
            end
        end

        if sol.converged
            prs[i] = sol.PR
            mdots[i] = sol.mdot
            powers[i] = sol.tau * omega
            converged[i] = true
            seed = (mdot=sol.mdot, ht_out=sol.ht_out, tau=sol.tau)
        else
            @warn "Operating-point solve failed at omega=$omega: retcode=$(sol.retcode)"
            seed = direct_try.branch_seed
        end
    end

    return (
        omegas=omegas,
        prs=prs,
        mdots=mdots,
        powers=powers,
        converged=converged,
        backoff_used=backoff_used,
    )
end

function _plot_compressor_operating_sweep_data(
    data;
    output_path::AbstractString="compressor_operating_sweep_via_operating_point_solver.png",
)
    p1 = plot(
        data.omegas,
        data.prs;
        xlabel="shaft speed (normalized)",
        ylabel="compression ratio Pt_out/Pt_in",
        title="Compression Ratio vs Shaft Speed",
        lw=2,
        marker=:circle,
        label=false,
    )

    p2 = plot(
        data.omegas,
        data.mdots;
        xlabel="shaft speed (normalized)",
        ylabel="mass flow rate mdot (kg/s)",
        title="Mass Flow vs Shaft Speed",
        lw=2,
        marker=:circle,
        label=false,
    )

    p3 = plot(
        data.omegas,
        data.powers;
        xlabel="shaft speed (normalized)",
        ylabel="power consumption (W)",
        title="Power Consumption vs Shaft Speed",
        lw=2,
        marker=:circle,
        label=false,
    )

    fig = plot(p1, p2, p3; layout=(3, 1), size=(900, 1100))
    savefig(fig, output_path)
    n_converged = count(data.converged)
    println("Converged points: $n_converged / $(length(data.omegas))")
    println("Saved operating-point sweep plot to: $output_path")
end

function plot_compressor_operating_sweep(
    map::TM.AbstractCompressorPerformanceMap;
    output_path::String="compressor_operating_sweep_via_operating_point_solver.png",
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
    target_pr::Float64=2.0,
    pr_backoff::Bool=true,
    backoff_min_pt_out::Union{Nothing,Float64}=nothing,
    backoff_max_pt_out::Union{Nothing,Float64}=nothing,
    backoff_pt_out_tol::Float64=50.0,
    backoff_max_iters::Int=24,
)
    data = sweep_compressor_operating_points(
        map;
        omega_min=omega_min,
        omega_max=omega_max,
        n_points=n_points,
        pt_in=pt_in,
        Tt_in=Tt_in,
        target_pr=target_pr,
        pr_backoff=pr_backoff,
        backoff_min_pt_out=backoff_min_pt_out,
        backoff_max_pt_out=backoff_max_pt_out,
        backoff_pt_out_tol=backoff_pt_out_tol,
        backoff_max_iters=backoff_max_iters,
    )
    _plot_compressor_operating_sweep_data(data; output_path=output_path)
end

function write_compressor_operating_sweep_csv(
    data;
    output_path::AbstractString="compressor_operating_sweep.csv",
)
    open(output_path, "w") do io
        println(io, "omega,PR,mdot,power,converged,backoff_used")
        for i in eachindex(data.omegas)
            println(
                io,
                "$(data.omegas[i]),$(data.prs[i]),$(data.mdots[i]),$(data.powers[i]),$(data.converged[i]),$(data.backoff_used[i])",
            )
        end
    end
    n_converged = count(data.converged)
    println("Converged points: $n_converged / $(length(data.omegas))")
    println("Saved operating-point sweep CSV to: $output_path")
    return output_path
end

function _main()
    settings = ArgParseSettings(
        prog="plot_compressor_operating_sweep.jl",
        description="Solve and plot compressor operating sweep from a compressor performance map file.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input compressor performance map (.toml)"
            required = true
        "--output"
            help = "output plot path"
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
        "--disable-pr-backoff"
            help = "disable pt_out backoff when target PR solve fails"
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
            help = "minimum normalized shaft speed"
            arg_type = Float64
            default = 0.6
        "--omega-max"
            help = "maximum normalized shaft speed"
            arg_type = Float64
            default = 1.0
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
    end

    parsed = parse_args(ARGS, settings)
    csv_mode = something(_parsed_opt(parsed, "csv", "csv"), false)
    map_group = _parsed_opt(parsed, "map_group", "map-group")
    output_default = csv_mode ? "compressor_operating_sweep.csv" : "compressor_operating_sweep.png"
    output = something(_parsed_opt(parsed, "output", "output"), output_default)
    target_pr = something(_parsed_opt(parsed, "target_pr", "target-pr"), 2.0)
    pr_backoff = !something(_parsed_opt(parsed, "disable_pr_backoff", "disable-pr-backoff"), false)
    backoff_min_pt_out = _parsed_opt(parsed, "backoff_min_pt_out", "backoff-min-pt-out")
    backoff_max_pt_out = _parsed_opt(parsed, "backoff_max_pt_out", "backoff-max-pt-out")
    backoff_pt_out_tol = something(_parsed_opt(parsed, "backoff_pt_out_tol", "backoff-pt-out-tol"), 50.0)
    backoff_max_iters = something(_parsed_opt(parsed, "backoff_max_iters", "backoff-max-iters"), 24)
    omega_min = something(_parsed_opt(parsed, "omega_min", "omega-min"), 0.6)
    omega_max = something(_parsed_opt(parsed, "omega_max", "omega-max"), 1.0)
    n_points = something(_parsed_opt(parsed, "n_points", "n-points"), 25)
    pt_in = something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0)
    tt_in = something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15)
    map = _load_compressor_map(parsed["map_path"]; group=map_group)
    data = sweep_compressor_operating_points(
        map;
        omega_min=omega_min,
        omega_max=omega_max,
        n_points=n_points,
        pt_in=pt_in,
        Tt_in=tt_in,
        target_pr=target_pr,
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
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
