"""
Compressor operating-point and sweep solve helpers.
"""

using NonlinearSolve
using ....Utility: bracket_bisect_roots, feasibility_backoff

"""
Solve a simple compressor operating point with four boundary conditions.

Fixed boundary conditions:
- `pt_in`
- `ht_in`
- `pt_out`
- `omega`

Solved unknowns:
- `mdot`
- `ht_out`
- `tau`
"""
function solve_compressor_operating_point(
    compressor_map::AbstractCompressorPerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    omega::Real,
    mdot_guess::Real,
    ht_out_guess::Real,
    tau_guess::Real,
    abstol::Real=1e-10,
    reltol::Real=1e-8,
    maxiters::Int=100,
    scaled_residuals::Bool=true,
    pressure_scale::Union{Nothing,Real}=nothing,
    enthalpy_scale::Union{Nothing,Real}=nothing,
    power_scale::Union{Nothing,Real}=nothing,
)
    scale_overrides = compressor_residual_scales(
        pt_in,
        ht_in,
        pt_out,
        ht_out_guess,
        mdot_guess,
        omega,
        tau_guess;
        pressure_scale=pressure_scale,
        enthalpy_scale=enthalpy_scale,
        power_scale=power_scale,
    )

    function f!(R, u, _)
        mdot = u[1]
        ht_out = u[2]
        tau = u[3]
        if scaled_residuals
            r_pout, r_dh_eff, r_P = compressor_residuals_scaled(
                compressor_map,
                eos,
                pt_in,
                ht_in,
                pt_out,
                ht_out,
                mdot,
                omega,
                tau;
                pressure_scale=scale_overrides.pressure_scale,
                enthalpy_scale=scale_overrides.enthalpy_scale,
                power_scale=scale_overrides.power_scale,
            )
            R[1] = r_pout
            R[2] = r_dh_eff
            R[3] = r_P
        else
            R_pout, R_dh_eff, R_P = compressor_residuals(
                compressor_map,
                eos,
                pt_in,
                ht_in,
                pt_out,
                ht_out,
                mdot,
                omega,
                tau,
            )
            R[1] = R_pout
            R[2] = R_dh_eff
            R[3] = R_P
        end
        return nothing
    end

    u0 = Float64[mdot_guess, ht_out_guess, tau_guess]
    prob = NonlinearProblem(f!, u0)
    sol = solve(
        prob,
        NewtonRaphson(; autodiff=AutoForwardDiff());
        abstol=abstol,
        reltol=reltol,
        maxiters=maxiters,
    )

    mdot = sol.u[1]
    ht_out = sol.u[2]
    tau = sol.u[3]

    residuals = compressor_residuals(
        compressor_map,
        eos,
        pt_in,
        ht_in,
        pt_out,
        ht_out,
        mdot,
        omega,
        tau,
    )

    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    map_vals = compressor_performance_map_from_stagnation(
        compressor_map,
        omega,
        mdot,
        Tt_in,
        pt_in,
    )

    return (
        mdot=mdot,
        ht_out=ht_out,
        tau=tau,
        PR=map_vals.PR,
        eta=map_vals.eta,
        residuals=residuals,
        retcode=sol.retcode,
        converged=(string(sol.retcode) == "Success"),
        solution=sol,
    )
end

function _compressor_flow_bounds(
    map::AbstractCompressorPerformanceMap,
    omega::Float64,
    Tt_in::Float64,
)
    domain = performance_map_domain(map)
    omega_corr = corrected_speed(omega, Tt_in, map)
    flow_range = domain.mdot_corr_flow_range
    m_surge = flow_range.surge(omega_corr)
    m_choke = flow_range.choke(omega_corr)
    return (min(m_surge, m_choke), max(m_surge, m_choke), omega_corr)
end

"""
    compressor_pr_roots(map; omega, Tt_in, target_pr, n_scan=401, pr_tol=1e-8, prior_roots=[])

Find all corrected-flow roots at fixed speed where `PR(omega, mdot_corr) = target_pr`.
Returns roots sorted by `mdot_corr` as `(mdot_corr, PR, eta)` tuples.
"""
function compressor_pr_roots(
    map::AbstractCompressorPerformanceMap;
    omega::Real,
    Tt_in::Real,
    target_pr::Real,
    n_scan::Int=401,
    pr_tol::Float64=1e-8,
    prior_roots::AbstractVector{<:Real}=Float64[],
)
    n_scan >= 5 || error("n_scan must be >= 5")
    target_pr > 0.0 || error("target_pr must be > 0")
    omega_f = Float64(omega)
    Tt_in_f = Float64(Tt_in)
    target_pr_f = Float64(target_pr)

    m_lo, m_hi, omega_corr = _compressor_flow_bounds(map, omega_f, Tt_in_f)
    f = mdot_corr -> begin
        vals = compressor_performance_map(map, omega_corr, mdot_corr)
        return vals.PR - target_pr_f
    end
    roots_corr = bracket_bisect_roots(
        f,
        (m_lo, m_hi);
        n_scan=n_scan,
        root_tol=pr_tol,
        prior_roots=prior_roots,
    )

    roots = NamedTuple[]
    for mdot_corr in roots_corr
        vals = compressor_performance_map(map, omega_corr, mdot_corr)
        push!(roots, (mdot_corr=mdot_corr, PR=vals.PR, eta=vals.eta))
    end
    return roots
end

function _compressor_pr_roots_with_backoff(
    map::AbstractCompressorPerformanceMap,
    omega::Float64,
    Tt_in::Float64,
    pt_in::Float64,
    target_pr::Float64;
    enabled::Bool=true,
    min_pt_out::Float64=pt_in,
    max_pt_out::Float64=pt_in * target_pr,
    pt_out_tol::Float64=50.0,
    max_iters::Int=24,
    prior_roots::AbstractVector{<:Real}=Float64[],
)
    target_pt_out = pt_in * target_pr
    lo_pt = max(min_pt_out, pt_in)
    hi_pt = min(max_pt_out, target_pt_out)
    lo_pt <= hi_pt || return (roots=NamedTuple[], PR=NaN, backoff_used=false, converged=false)

    evaluate_at_pt = pt -> begin
        pr = pt / pt_in
        compressor_pr_roots(
            map;
            omega=omega,
            Tt_in=Tt_in,
            target_pr=pr,
            prior_roots=prior_roots,
        )
    end

    backoff = feasibility_backoff(
        evaluate_at_pt,
        target_pt_out;
        enabled=enabled,
        min_value=lo_pt,
        max_value=hi_pt,
        value_tol=pt_out_tol,
        max_iters=max_iters,
        n_probe=33,
        is_feasible=(roots -> !isempty(roots)),
    )

    backoff.converged || return (roots=NamedTuple[], PR=NaN, backoff_used=false, converged=false)
    return (
        roots=backoff.result,
        PR=(backoff.value / pt_in),
        backoff_used=backoff.used_backoff,
        converged=true,
    )
end

function _compressor_select_roots(roots::Vector{NamedTuple}, branch::Symbol)
    isempty(roots) && return NamedTuple[]
    if branch == :low
        return [first(roots)]
    elseif branch == :high
        return [last(roots)]
    else
        return roots
    end
end

function _compressor_root_outputs(
    root::NamedTuple,
    eos::Fluids.AbstractEOS,
    map::AbstractCompressorPerformanceMap,
    pt_in::Float64,
    ht_in::Float64,
    Tt_in::Float64,
    omega::Float64,
)
    corr_to_phys_scale = (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
    mdot = root.mdot_corr * corr_to_phys_scale
    pt_out = pt_in * root.PR
    eta_safe = max(root.eta, 1e-6)
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    ht_out = ht_in + (h2s - ht_in) / eta_safe
    tau = mdot * (ht_out - ht_in) / max(omega, 1e-12)
    power = tau * omega
    return (PR=root.PR, eta=root.eta, mdot=mdot, ht_out=ht_out, tau=tau, power=power)
end

"""
    solve_compressor_operating_sweep(map, eos; ...)

Sweep compressor operating points across `omega` for fixed `pt_in`, `Tt_in`, and target `PR`.
Supports branch selection:
- `:high` (default): highest-flow root
- `:low`: lowest-flow root
- `:all`: all roots at each speed (CSV-oriented output)

Returns either:
- `mode=:single` with one row per speed (`:low`/`:high`)
- `mode=:all` with `rows` containing all roots
"""
function solve_compressor_operating_sweep(
    map::AbstractCompressorPerformanceMap,
    eos::Fluids.AbstractEOS;
    omega_min::Real=0.6,
    omega_max::Real=1.0,
    n_points::Int=25,
    pt_in::Real=101_325.0,
    Tt_in::Real=288.15,
    target_pr::Real=2.0,
    branch::Symbol=:high,
    pr_backoff::Bool=true,
    backoff_min_pt_out::Union{Nothing,Real}=nothing,
    backoff_max_pt_out::Union{Nothing,Real}=nothing,
    backoff_pt_out_tol::Real=50.0,
    backoff_max_iters::Int=24,
)
    target_pr_f = Float64(target_pr)
    target_pr_f > 1.0 || error("target_pr must be > 1.0")
    branch in (:low, :high, :all) || error("branch must be one of: low|high|all")

    omega_min_f = Float64(omega_min)
    omega_max_f = Float64(omega_max)
    pt_in_f = Float64(pt_in)
    Tt_in_f = Float64(Tt_in)
    ht_in = Fluids.enthalpy_from_temperature(eos, Tt_in_f)
    omegas = collect(range(omega_min_f, omega_max_f, length=n_points))

    min_pt_out = isnothing(backoff_min_pt_out) ? pt_in_f : Float64(backoff_min_pt_out)
    max_pt_out = isnothing(backoff_max_pt_out) ? (pt_in_f * target_pr_f) : Float64(backoff_max_pt_out)

    if branch == :all
        rows = NamedTuple[]
        prior_roots = Float64[]
        for omega in omegas
            found = _compressor_pr_roots_with_backoff(
                map,
                omega,
                Tt_in_f,
                pt_in_f,
                target_pr_f;
                enabled=pr_backoff,
                min_pt_out=min_pt_out,
                max_pt_out=max_pt_out,
                pt_out_tol=Float64(backoff_pt_out_tol),
                max_iters=backoff_max_iters,
                prior_roots=prior_roots,
            )

            if !found.converged || isempty(found.roots)
                push!(rows, (
                    omega=omega,
                    branch_id=0,
                    PR=NaN,
                    eta=NaN,
                    mdot=NaN,
                    power=NaN,
                    converged=false,
                    backoff_used=false,
                ))
                continue
            end

            selected = _compressor_select_roots(found.roots, :all)
            for (k, root) in enumerate(selected)
                out = _compressor_root_outputs(root, eos, map, pt_in_f, ht_in, Tt_in_f, omega)
                push!(rows, (
                    omega=omega,
                    branch_id=k,
                    PR=out.PR,
                    eta=out.eta,
                    mdot=out.mdot,
                    power=out.power,
                    converged=true,
                    backoff_used=found.backoff_used,
                ))
            end
            prior_roots = [r.mdot_corr for r in found.roots]
        end
        return (mode=:all, branch=:all, rows=rows)
    end

    prs = fill(NaN, n_points)
    etas = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    converged = fill(false, n_points)
    backoff_used = fill(false, n_points)
    prior_roots = Float64[]

    for (i, omega) in enumerate(omegas)
        found = _compressor_pr_roots_with_backoff(
            map,
            omega,
            Tt_in_f,
            pt_in_f,
            target_pr_f;
            enabled=pr_backoff,
            min_pt_out=min_pt_out,
            max_pt_out=max_pt_out,
            pt_out_tol=Float64(backoff_pt_out_tol),
            max_iters=backoff_max_iters,
            prior_roots=prior_roots,
        )

        if !found.converged || isempty(found.roots)
            continue
        end

        root = only(_compressor_select_roots(found.roots, branch))
        out = _compressor_root_outputs(root, eos, map, pt_in_f, ht_in, Tt_in_f, omega)
        prs[i] = out.PR
        etas[i] = out.eta
        mdots[i] = out.mdot
        powers[i] = out.power
        converged[i] = true
        backoff_used[i] = found.backoff_used
        prior_roots = [r.mdot_corr for r in found.roots]
    end

    return (
        mode=:single,
        branch=branch,
        omegas=omegas,
        prs=prs,
        etas=etas,
        mdots=mdots,
        powers=powers,
        converged=converged,
        backoff_used=backoff_used,
    )
end

function solve_compressor_operating_sweep(
    map::AbstractCompressorPerformanceMap;
    kwargs...,
)
    eos = Fluids.ideal_EOS()[:air]
    return solve_compressor_operating_sweep(map, eos; kwargs...)
end
