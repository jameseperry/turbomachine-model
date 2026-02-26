#!/usr/bin/env julia

using TurboMachineModel
using Plots

const Fluids = TurboMachineModel.Physics.Fluids
const TT = TurboMachineModel.Physics.Turbomachine.Turbine

function _state_from_pt_out(
    map::TT.TabulatedTurbinePerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    pt_out::Real,
)
    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    map_vals = TT.turbine_performance_map_from_stagnation(map, omega, pt_in, pt_out, Tt_in)
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    ht_out = ht_in - map_vals.eta * (ht_in - h2s)
    power = map_vals.mdot * (ht_out - ht_in)
    return (
        PR_turb=map_vals.PR_turb,
        eta=map_vals.eta,
        mdot=map_vals.mdot,
        ht_out=ht_out,
        power=power,
    )
end

function _solve_pt_out_at_speed(
    map::TT.TabulatedTurbinePerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    tau_load::Real,
    pt_out_guess::Real,
)
    function mismatch(pt_out::Float64)
        state = _state_from_pt_out(map, eos, pt_in, ht_in, omega, pt_out)
        return tau_load * omega - state.power
    end

    pr_min = first(map.pr_turb_grid)
    pr_max = last(map.pr_turb_grid)
    pt_out_min = pt_in / pr_max
    pt_out_max = pt_in / pr_min
    tau_ref = max(abs(tau_load * omega), 1.0)
    mismatch_tol = max(1e-8 * tau_ref, 1e-4)
    pt_tol = 1e-10 * max(pt_in, 1.0)

    scan_pts = collect(range(pt_out_min, pt_out_max, length=31))
    if pt_out_min <= pt_out_guess <= pt_out_max
        push!(scan_pts, pt_out_guess)
        sort!(scan_pts)
    end

    a = NaN
    b = NaN
    fa = NaN
    fb = NaN
    best_pt = NaN
    best_abs = Inf

    for i in eachindex(scan_pts)
        pt = scan_pts[i]
        f = mismatch(pt)
        if !isfinite(f)
            continue
        end
        abs_f = abs(f)
        if abs_f < best_abs
            best_abs = abs_f
            best_pt = pt
        end
        if i < length(scan_pts)
            pt2 = scan_pts[i + 1]
            f2 = mismatch(pt2)
            if isfinite(f2) && signbit(f) != signbit(f2)
                a = pt
                b = pt2
                fa = f
                fb = f2
                break
            end
        end
    end

    if isfinite(best_abs) && best_abs <= mismatch_tol
        return (converged=true, pt_out=best_pt, reason="")
    end

    if !(isfinite(a) && isfinite(b))
        return (converged=false, pt_out=best_pt, reason="no sign-change bracket in map range")
    end

    for _ in 1:80
        mid = 0.5 * (a + b)
        fm = mismatch(mid)
        if !isfinite(fm)
            return (converged=false, pt_out=best_pt, reason="non-finite mismatch during bisection")
        end
        if abs(fm) <= mismatch_tol || abs(b - a) <= pt_tol
            return (converged=true, pt_out=mid, reason="")
        end
        if signbit(fm) == signbit(fa)
            a = mid
            fa = fm
        else
            b = mid
            fb = fm
        end
    end

    return (converged=false, pt_out=best_pt, reason="bisection max iterations")
end

function sweep_turbine_operating_points(;
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    tau_load::Float64=-5.0e5,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
)
    map = TT.demo_turbine_performance_map()
    eos = Fluids.ideal_EOS()[:air]
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in)

    omegas = collect(range(omega_min, omega_max, length=n_points))
    prs = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    converged = fill(false, n_points)

    pr_mid = 0.5 * (first(map.pr_turb_grid) + last(map.pr_turb_grid))
    pt_out_guess = pt_in / pr_mid

    for (i, omega) in enumerate(omegas)
        sol = _solve_pt_out_at_speed(map, eos, pt_in, ht_in, omega, tau_load, pt_out_guess)
        if sol.converged
            pt_out = sol.pt_out
            state = _state_from_pt_out(map, eos, pt_in, ht_in, omega, pt_out)
            prs[i] = state.PR_turb
            mdots[i] = state.mdot
            powers[i] = -state.power
            converged[i] = true
            pt_out_guess = pt_out
        else
            @warn "Operating-point solve failed at omega=$omega: $(sol.reason)"
        end
    end

    return (omegas=omegas, prs=prs, mdots=mdots, powers=powers, converged=converged)
end

function plot_turbine_operating_sweep(;
    output_path::String="turbine_operating_sweep.png",
    tau_load::Float64=-5.0e5,
)
    data = sweep_turbine_operating_points(tau_load=tau_load)

    p1 = plot(
        data.omegas,
        data.prs;
        xlabel="shaft speed (normalized)",
        ylabel="expansion ratio Pt_in/Pt_out",
        title="Expansion Ratio vs Shaft Speed",
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
        ylabel="shaft power output (-tau*omega) (W)",
        title="Power Output vs Shaft Speed",
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
    output_path = length(ARGS) >= 1 ? ARGS[1] : "turbine_operating_sweep.png"
    tau_load = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : -5.0e5
    plot_turbine_operating_sweep(output_path=output_path, tau_load=tau_load)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
