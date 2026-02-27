#!/usr/bin/env julia

using ArgParse
using TurboMachineModel
using Plots

const Fluids = TurboMachineModel.Physics.Fluids
const TM = TurboMachineModel.Physics.Turbomachine.Compressor
const U = TurboMachineModel.Utility

function _infer_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".toml"
        return :toml
    end
    error("unsupported map extension $(ext) for path $(path); expected .toml")
end

function _load_compressor_map(path::AbstractString; group::AbstractString="compressor_map")
    _infer_format(path)
    return TM.read_toml(TM.TabulatedCompressorPerformanceMap, path; group=group)
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _inner_operating_point(
    map::TM.TabulatedCompressorPerformanceMap,
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

function _tau_mismatch(
    map::TM.TabulatedCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    tau_load::Float64,
    pt_out::Float64,
    guess,
)
    sol = _inner_operating_point(
        map,
        eos,
        pt_in,
        ht_in,
        omega,
        pt_out,
        guess.mdot,
        guess.ht_out,
        guess.tau,
    )
    if !sol.converged
        return (ok=false, mismatch=NaN, sol=sol)
    end
    return (ok=true, mismatch=sol.tau - tau_load, sol=sol)
end

function _solve_at_speed_via_operating_point(
    map::TM.TabulatedCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    tau_load::Float64,
    seed,
)
    pr_min = max(1.01, minimum(U.table_values(map.pr_map)) * 0.9)
    pr_max = maximum(U.table_values(map.pr_map)) * 1.1
    pt_min = pt_in * pr_min
    pt_max = pt_in * pr_max

    n_scan = 17
    scan_pts = collect(range(pt_min, pt_max, length=n_scan))
    scan_vals = Vector{Tuple{Float64,Float64,Any}}()
    local_guess = (mdot=seed.mdot, ht_out=seed.ht_out, tau=seed.tau)

    for pt_out in scan_pts
        eval = _tau_mismatch(map, eos, pt_in, ht_in, omega, tau_load, pt_out, local_guess)
        if eval.ok
            push!(scan_vals, (pt_out, eval.mismatch, eval.sol))
            local_guess = (mdot=eval.sol.mdot, ht_out=eval.sol.ht_out, tau=eval.sol.tau)
        end
    end

    length(scan_vals) >= 2 || return (converged=false, reason="no converged inner solves")

    bracket_idx = nothing
    for k in 1:(length(scan_vals) - 1)
        left = scan_vals[k]
        right = scan_vals[k + 1]
        sign(left[2]) == sign(right[2]) && continue
        bracket_idx = k
        break
    end
    bracket_idx === nothing && return (converged=false, reason="no bracket found")

    a_pt, a_mis, a_sol = scan_vals[bracket_idx]
    b_pt, b_mis, b_sol = scan_vals[bracket_idx + 1]
    tau_tol = max(1e-8 * max(abs(tau_load), 1.0), 1e-4)
    pt_tol = 1e-8 * max(pt_in, 1.0)

    abs(a_mis) < abs(b_mis) ? (best_pt, best_mis, best_sol) = (a_pt, a_mis, a_sol) :
                              (best_pt, best_mis, best_sol) = (b_pt, b_mis, b_sol)

    for _ in 1:40
        mid_pt = 0.5 * (a_pt + b_pt)
        mid_guess = abs(a_mis) <= abs(b_mis) ?
                    (mdot=a_sol.mdot, ht_out=a_sol.ht_out, tau=a_sol.tau) :
                    (mdot=b_sol.mdot, ht_out=b_sol.ht_out, tau=b_sol.tau)

        mid_eval = _tau_mismatch(map, eos, pt_in, ht_in, omega, tau_load, mid_pt, mid_guess)
        if !mid_eval.ok
            alt_guess = abs(a_mis) > abs(b_mis) ?
                        (mdot=a_sol.mdot, ht_out=a_sol.ht_out, tau=a_sol.tau) :
                        (mdot=b_sol.mdot, ht_out=b_sol.ht_out, tau=b_sol.tau)
            mid_eval = _tau_mismatch(map, eos, pt_in, ht_in, omega, tau_load, mid_pt, alt_guess)
            if !mid_eval.ok
                continue
            end
        end

        mid_mis = mid_eval.mismatch
        mid_sol = mid_eval.sol
        if abs(mid_mis) < abs(best_mis)
            best_pt, best_mis, best_sol = mid_pt, mid_mis, mid_sol
        end

        if abs(mid_mis) <= tau_tol || abs(b_pt - a_pt) <= pt_tol
            return (
                converged=true,
                pt_out=mid_pt,
                mdot=mid_sol.mdot,
                ht_out=mid_sol.ht_out,
                tau=mid_sol.tau,
                PR=mid_pt / pt_in,
            )
        end

        if sign(mid_mis) == sign(a_mis)
            a_pt, a_mis, a_sol = mid_pt, mid_mis, mid_sol
        else
            b_pt, b_mis, b_sol = mid_pt, mid_mis, mid_sol
        end
    end

    if abs(best_mis) <= tau_tol
        return (
            converged=true,
            pt_out=best_pt,
            mdot=best_sol.mdot,
            ht_out=best_sol.ht_out,
            tau=best_sol.tau,
            PR=best_pt / pt_in,
        )
    end

    return (converged=false, reason="bisection did not converge")
end

function sweep_compressor_operating_points(
    map::TM.TabulatedCompressorPerformanceMap;
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    tau_load::Float64=1.2e6,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
)
    eos = Fluids.ideal_EOS()[:air]
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)

    omegas = collect(range(omega_min, omega_max, length=n_points))
    prs = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    converged = fill(false, n_points)

    mdot_scale = (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
    mdot0 = 0.5 * (first(U.table_ygrid(map.pr_map)) + last(U.table_ygrid(map.pr_map))) * mdot_scale
    map0 = TM.compressor_performance_map_from_stagnation(map, omegas[1], mdot0, Tt_in, pt_in)
    pt_out0 = pt_in * map0.PR
    h2s0 = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out0)
    ht_out0 = ht_in + (h2s0 - ht_in) / map0.eta
    tau0 = mdot0 * (ht_out0 - ht_in) / omegas[1]
    seed = (mdot=mdot0, ht_out=ht_out0, tau=tau0)

    for (i, omega) in enumerate(omegas)
        sol = _solve_at_speed_via_operating_point(
            map,
            eos,
            pt_in,
            ht_in,
            omega,
            tau_load,
            seed,
        )
        if sol.converged
            prs[i] = sol.PR
            mdots[i] = sol.mdot
            powers[i] = sol.tau * omega
            converged[i] = true
            seed = (mdot=sol.mdot, ht_out=sol.ht_out, tau=sol.tau)
        else
            @warn "Operating-point solve failed at omega=$omega: $(sol.reason)"
        end
    end

    return (omegas=omegas, prs=prs, mdots=mdots, powers=powers, converged=converged)
end

function plot_compressor_operating_sweep(
    map::TM.TabulatedCompressorPerformanceMap;
    output_path::String="compressor_operating_sweep_via_operating_point_solver.png",
    tau_load::Float64=1.2e6,
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
)
    data = sweep_compressor_operating_points(
        map;
        tau_load=tau_load,
        omega_min=omega_min,
        omega_max=omega_max,
        n_points=n_points,
        pt_in=pt_in,
        Tt_in=Tt_in,
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

function _main()
    settings = ArgParseSettings(
        prog="plot_compressor_operating_sweep.jl",
        description="Solve and plot compressor operating sweep from a tabulated performance map file.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input compressor performance map (.toml)"
            required = true
        "--output"
            help = "output plot path"
            arg_type = String
            default = "compressor_operating_sweep.png"
        "--map-group"
            help = "input map group/table"
            arg_type = String
            default = "compressor_map"
        "--tau-load"
            help = "load torque (N*m)"
            arg_type = Float64
            default = 1.2e6
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
    map_group = something(_parsed_opt(parsed, "map_group", "map-group"), "compressor_map")
    output = something(_parsed_opt(parsed, "output", "output"), "compressor_operating_sweep.png")
    tau_load = something(_parsed_opt(parsed, "tau_load", "tau-load"), 1.2e6)
    omega_min = something(_parsed_opt(parsed, "omega_min", "omega-min"), 0.6)
    omega_max = something(_parsed_opt(parsed, "omega_max", "omega-max"), 1.0)
    n_points = something(_parsed_opt(parsed, "n_points", "n-points"), 25)
    pt_in = something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0)
    tt_in = something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15)
    map = _load_compressor_map(parsed["map_path"]; group=map_group)
    plot_compressor_operating_sweep(
        map;
        output_path=output,
        tau_load=tau_load,
        omega_min=omega_min,
        omega_max=omega_max,
        n_points=n_points,
        pt_in=pt_in,
        Tt_in=tt_in,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
