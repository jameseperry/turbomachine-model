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
    elseif ext == ".h5" || ext == ".hdf5"
        return :hdf5
    end
    error("unsupported map extension $(ext) for path $(path); expected .toml/.h5/.hdf5")
end

function _load_turbine_map(path::AbstractString; group::AbstractString="turbine_map")
    format = _infer_format(path)
    if format == :toml
        return TT.read_toml(TT.TabulatedTurbinePerformanceMap, path; group=group)
    end
    return TT.read_hdf5(TT.TabulatedTurbinePerformanceMap, path; group=group)
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

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

    pr_min = first(U.table_ygrid(map.mdot_corr_map))
    pr_max = last(U.table_ygrid(map.mdot_corr_map))
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

function sweep_turbine_operating_points(
    map::TT.TabulatedTurbinePerformanceMap;
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    tau_load::Float64=-5.0e5,
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

    pr_mid = 0.5 * (first(U.table_ygrid(map.mdot_corr_map)) + last(U.table_ygrid(map.mdot_corr_map)))
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

function plot_turbine_operating_sweep(
    map::TT.TabulatedTurbinePerformanceMap;
    output_path::String="turbine_operating_sweep.png",
    tau_load::Float64=-5.0e5,
    omega_min::Float64=0.6,
    omega_max::Float64=1.0,
    n_points::Int=25,
    pt_in::Float64=101_325.0,
    Tt_in::Float64=288.15,
)
    data = sweep_turbine_operating_points(
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
    settings = ArgParseSettings(
        prog="plot_turbine_operating_sweep.jl",
        description="Solve and plot turbine operating sweep from a tabulated performance map file.",
    )
    @add_arg_table! settings begin
        "map_path"
            help = "input turbine performance map (.toml/.h5/.hdf5)"
            required = true
        "--output"
            help = "output plot path"
            arg_type = String
            default = "turbine_operating_sweep.png"
        "--map-group"
            help = "input map group/table"
            arg_type = String
            default = "turbine_map"
        "--tau-load"
            help = "load torque (N*m), typically negative for power extraction"
            arg_type = Float64
            default = -5.0e5
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
    map_group = something(_parsed_opt(parsed, "map_group", "map-group"), "turbine_map")
    output = something(_parsed_opt(parsed, "output", "output"), "turbine_operating_sweep.png")
    tau_load = something(_parsed_opt(parsed, "tau_load", "tau-load"), -5.0e5)
    omega_min = something(_parsed_opt(parsed, "omega_min", "omega-min"), 0.6)
    omega_max = something(_parsed_opt(parsed, "omega_max", "omega-max"), 1.0)
    n_points = something(_parsed_opt(parsed, "n_points", "n-points"), 25)
    pt_in = something(_parsed_opt(parsed, "pt_in", "pt-in"), 101_325.0)
    tt_in = something(_parsed_opt(parsed, "tt_in", "tt-in"), 288.15)
    map = _load_turbine_map(parsed["map_path"]; group=map_group)
    plot_turbine_operating_sweep(
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
