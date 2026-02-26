#!/usr/bin/env julia

using TurboMachineModel
using NonlinearSolve
using Plots

const Fluids = TurboMachineModel.Physics.Fluids
const TM = TurboMachineModel.Physics.Turbomachine.Compressor

function _state_from_mdot(
    map::TM.TabulatedCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    mdot::Real,
)
    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    map_vals = TM.compressor_performance_map_from_stagnation(map, omega, mdot, Tt_in, pt_in)
    pt_out = map_vals.PR * pt_in
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    ht_out = ht_in + (h2s - ht_in) / map_vals.eta
    power = mdot * (ht_out - ht_in)
    return (
        PR=map_vals.PR,
        eta=map_vals.eta,
        pt_out=pt_out,
        ht_out=ht_out,
        power=power,
    )
end

function _solve_mdot_at_speed(
    map::TM.TabulatedCompressorPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    omega::Real,
    tau_load::Real,
    mdot_guess::Real,
)
    function residual(u, _)
        mdot = u[1]
        state = _state_from_mdot(map, eos, pt_in, ht_in, omega, mdot)
        return [tau_load * omega - state.power]
    end

    function solve_with_guess(guess::Float64)
        prob = NonlinearProblem(residual, [guess])
        return solve(
            prob,
            NewtonRaphson(; autodiff=AutoForwardDiff());
            abstol=1e-10,
            reltol=1e-8,
            maxiters=100,
        )
    end

    sol = solve_with_guess(mdot_guess)
    string(sol.retcode) == "Success" && return sol

    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    mdot_scale = (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
    corrected_min = first(map.mdot_corr_grid)
    corrected_max = last(map.mdot_corr_grid)
    corrected_mid = 0.5 * (corrected_min + corrected_max)
    fallback_guesses = (
        corrected_mid * mdot_scale,
        corrected_min * mdot_scale,
        corrected_max * mdot_scale,
        0.8 * corrected_mid * mdot_scale,
        1.2 * corrected_mid * mdot_scale,
    )

    for guess in fallback_guesses
        sol_try = solve_with_guess(guess)
        if string(sol_try.retcode) == "Success"
            return sol_try
        end
    end

    return sol
end

function sweep_compressor_operating_points(;
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    tau_load::Float64=3.0e5,
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
    mdot_guess = 0.5 * (first(map.mdot_corr_grid) + last(map.mdot_corr_grid)) * mdot_scale

    for (i, omega) in enumerate(omegas)
        sol = _solve_mdot_at_speed(map, eos, pt_in, ht_in, omega, tau_load, mdot_guess)
        if string(sol.retcode) == "Success"
            mdot = sol.u[1]
            state = _state_from_mdot(map, eos, pt_in, ht_in, omega, mdot)
            prs[i] = state.PR
            mdots[i] = mdot
            powers[i] = state.power
            converged[i] = true
            mdot_guess = mdot
        else
            @warn "Operating-point solve failed at omega=$omega with retcode=$(sol.retcode)"
        end
    end

    return (omegas=omegas, prs=prs, mdots=mdots, powers=powers, converged=converged)
end

function plot_compressor_operating_sweep(;
    output_path::String="compressor_operating_sweep.png",
    tau_load::Float64=3.0e5,
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
        ylabel="power consumption tau*omega (W)",
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
    output_path = length(ARGS) >= 1 ? ARGS[1] : "compressor_operating_sweep.png"
    tau_load = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 3.0e5
    plot_compressor_operating_sweep(output_path=output_path, tau_load=tau_load)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end
