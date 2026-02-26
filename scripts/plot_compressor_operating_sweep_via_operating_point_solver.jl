#!/usr/bin/env julia

using TurboMachineModel
using Plots

const Fluids = TurboMachineModel.Physics.Fluids
const TM = TurboMachineModel.Physics.Turbomachine.Compressor

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
    pr_min = max(1.01, minimum(map.pr_table) * 0.9)
    pr_max = maximum(map.pr_table) * 1.1
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

function sweep_compressor_operating_points(;
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    tau_load::Float64=1.2e6,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
)
    map = TM.demo_compressor_performance_map()
    eos = Fluids.ideal_EOS()[:air]
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)

    omegas = collect(range(omega_min, omega_max, length=n_points))
    prs = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    converged = fill(false, n_points)

    mdot_scale = (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
    mdot0 = 0.5 * (first(map.mdot_corr_grid) + last(map.mdot_corr_grid)) * mdot_scale
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

function plot_compressor_operating_sweep(;
    output_path::String="compressor_operating_sweep_via_operating_point_solver.png",
    tau_load::Float64=1.2e6,
)
    data = sweep_compressor_operating_points(tau_load=tau_load)

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
    output_path = length(ARGS) >= 1 ? ARGS[1] :
                  "compressor_operating_sweep_via_operating_point_solver.png"
    tau_load = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.2e6
    plot_compressor_operating_sweep(output_path=output_path, tau_load=tau_load)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
