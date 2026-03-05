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
    Pt_in::Float64,
)
    domain = performance_map_domain(map, Tt_in, Pt_in)
    flow_range = domain.mdot_flow_range
    m_surge = flow_range.surge(omega)
    m_choke = flow_range.choke(omega)
    return (min(m_surge, m_choke), max(m_surge, m_choke))
end

function _compressor_condition_diagnostics(
    map::AbstractCompressorPerformanceMap,
    omega::Float64,
    Tt_in::Float64,
    pt_in::Float64,
    target_pr::Float64;
    n_scan::Int=401,
)
    domain = performance_map_domain(map, Tt_in, pt_in)
    flow_range = domain.mdot_flow_range
    mdot_surge = flow_range.surge(omega)
    mdot_choke = flow_range.choke(omega)
    flow_mdot_min = min(mdot_surge, mdot_choke)
    flow_mdot_max = max(mdot_surge, mdot_choke)

    pr_at = mdot -> begin
        vals = compressor_performance_map_from_stagnation(map, omega, mdot, Tt_in, pt_in)
        return vals.PR
    end

    pr_surge = pr_at(mdot_surge)
    pr_choke = pr_at(mdot_choke)

    pr_max = -Inf
    mdot_at_pr_max = NaN
    root_bracketed_by_scan = false
    has_prev = false
    f_prev = NaN

    n_scan_eff = max(n_scan, 5)
    for mdot in range(flow_mdot_min, flow_mdot_max, length=n_scan_eff)
        pr = pr_at(mdot)
        isfinite(pr) || continue

        if pr > pr_max
            pr_max = pr
            mdot_at_pr_max = mdot
        end

        f = pr - target_pr
        if has_prev && isfinite(f_prev)
            if f == 0.0 || f_prev == 0.0 || signbit(f) != signbit(f_prev)
                root_bracketed_by_scan = true
            end
        end
        f_prev = f
        has_prev = true
    end

    if !isfinite(pr_max)
        return (
            flow_mdot_min=flow_mdot_min,
            flow_mdot_max=flow_mdot_max,
            pr_surge=pr_surge,
            pr_choke=pr_choke,
            pr_max=NaN,
            mdot_at_pr_max=NaN,
            pr_peak_boundary="unknown",
            target_feasible_by_scan=false,
            root_bracketed_by_scan=false,
            infeasibility_hint="no_finite_pr_on_flow_interval",
        )
    end

    dist_to_surge = abs(mdot_at_pr_max - mdot_surge)
    dist_to_choke = abs(mdot_at_pr_max - mdot_choke)
    pr_peak_boundary = dist_to_surge <= dist_to_choke ? "surge" : "choke"
    target_feasible_by_scan = (pr_max >= target_pr)

    infeasibility_hint = if target_feasible_by_scan
        root_bracketed_by_scan ? "target_pr_bracketed" : "target_pr_not_bracketed_possible_tangent"
    else
        pr_peak_boundary == "surge" ? "target_pr_above_map_max_near_surge" : "target_pr_above_map_max_near_choke"
    end

    return (
        flow_mdot_min=flow_mdot_min,
        flow_mdot_max=flow_mdot_max,
        pr_surge=pr_surge,
        pr_choke=pr_choke,
        pr_max=pr_max,
        mdot_at_pr_max=mdot_at_pr_max,
        pr_peak_boundary=pr_peak_boundary,
        target_feasible_by_scan=target_feasible_by_scan,
        root_bracketed_by_scan=root_bracketed_by_scan,
        infeasibility_hint=infeasibility_hint,
    )
end

"""
    compressor_pr_roots(map; omega, Tt_in, Pt_in, target_pr, n_scan=401, pr_tol=1e-8, prior_roots=[])

Find all physical-flow roots at fixed speed where `PR(omega, mdot) = target_pr`.
Returns roots sorted by `mdot` as `(mdot, PR, eta)` tuples.
"""
function compressor_pr_roots(
    map::AbstractCompressorPerformanceMap;
    omega::Real,
    Tt_in::Real,
    Pt_in::Real,
    target_pr::Real,
    n_scan::Int=401,
    pr_tol::Float64=1e-8,
    prior_roots::AbstractVector{<:Real}=Float64[],
)
    n_scan >= 5 || error("n_scan must be >= 5")
    target_pr > 0.0 || error("target_pr must be > 0")
    omega_f = Float64(omega)
    Tt_in_f = Float64(Tt_in)
    Pt_in_f = Float64(Pt_in)
    target_pr_f = Float64(target_pr)

    m_lo, m_hi = _compressor_flow_bounds(map, omega_f, Tt_in_f, Pt_in_f)
    f = mdot -> begin
        vals = compressor_performance_map_from_stagnation(map, omega_f, mdot, Tt_in_f, Pt_in_f)
        return vals.PR - target_pr_f
    end
    roots_mdot = bracket_bisect_roots(
        f,
        (m_lo, m_hi);
        n_scan=n_scan,
        root_tol=pr_tol,
        prior_roots=prior_roots,
    )

    roots = NamedTuple[]
    for mdot in roots_mdot
        vals = compressor_performance_map_from_stagnation(map, omega_f, mdot, Tt_in_f, Pt_in_f)
        push!(roots, (mdot=mdot, PR=vals.PR, eta=vals.eta))
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
    lo_pt <= hi_pt || return (
        roots=NamedTuple[],
        PR=NaN,
        pt_out=NaN,
        backoff_used=false,
        converged=false,
        target_pr=target_pr,
        target_pt_out=target_pt_out,
        n_target_roots=0,
        failure_reason=:invalid_backoff_range,
    )

    evaluate_at_pt = pt -> begin
        pr = pt / pt_in
        compressor_pr_roots(
            map;
            omega=omega,
            Tt_in=Tt_in,
            Pt_in=pt_in,
            target_pr=pr,
            prior_roots=prior_roots,
        )
    end

    roots_hi = evaluate_at_pt(hi_pt)
    if !isempty(roots_hi)
        return (
            roots=roots_hi,
            PR=(hi_pt / pt_in),
            pt_out=hi_pt,
            backoff_used=(hi_pt < target_pt_out),
            converged=true,
            target_pr=target_pr,
            target_pt_out=target_pt_out,
            n_target_roots=length(roots_hi),
            failure_reason=:none,
        )
    end

    enabled || return (
        roots=NamedTuple[],
        PR=NaN,
        pt_out=NaN,
        backoff_used=false,
        converged=false,
        target_pr=target_pr,
        target_pt_out=target_pt_out,
        n_target_roots=length(roots_hi),
        failure_reason=:target_infeasible_backoff_disabled,
    )

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

    backoff.converged || return (
        roots=NamedTuple[],
        PR=NaN,
        pt_out=NaN,
        backoff_used=false,
        converged=false,
        target_pr=target_pr,
        target_pt_out=target_pt_out,
        n_target_roots=length(roots_hi),
        failure_reason=:no_feasible_backoff_solution,
    )
    return (
        roots=backoff.result,
        PR=(backoff.value / pt_in),
        pt_out=backoff.value,
        backoff_used=backoff.used_backoff,
        converged=true,
        target_pr=target_pr,
        target_pt_out=target_pt_out,
        n_target_roots=length(roots_hi),
        failure_reason=:none,
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

function _compressor_operating_margins(
    map::AbstractCompressorPerformanceMap,
    omega::Float64,
    mdot::Float64,
    Tt_in::Float64,
    pt_in::Float64,
)
    domain = performance_map_domain(map, Tt_in, pt_in)
    flow_range = domain.mdot_flow_range
    mdot_surge = flow_range.surge(omega)
    mdot_choke = flow_range.choke(omega)
    span = max(abs(mdot_choke - mdot_surge), 1e-12)
    surge_margin = mdot - mdot_surge
    choke_margin = mdot_choke - mdot

    map_vals = compressor_performance_map_from_stagnation(map, omega, mdot, Tt_in, pt_in)
    stall_flag = hasproperty(map_vals, :stall) ? Bool(map_vals.stall) : false
    choke_flag = hasproperty(map_vals, :choke) ? Bool(map_vals.choke) : false

    return (
        mdot_surge=mdot_surge,
        mdot_choke=mdot_choke,
        surge_margin=surge_margin,
        choke_margin=choke_margin,
        surge_margin_norm=surge_margin / span,
        choke_margin_norm=choke_margin / span,
        stall_flag=stall_flag,
        choke_flag=choke_flag,
    )
end

function _compressor_root_outputs(
    root::NamedTuple,
    eos::Fluids.AbstractEOS,
    map::AbstractCompressorPerformanceMap,
    pt_in::Float64,
    ht_in::Float64,
    omega::Float64,
    Tt_in::Float64,
)
    mdot = root.mdot
    pt_out = pt_in * root.PR
    eta_safe = max(root.eta, 1e-6)
    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
    ht_out = ht_in + (h2s - ht_in) / eta_safe
    tau = mdot * (ht_out - ht_in) / max(omega, 1e-12)
    power = tau * omega
    margins = _compressor_operating_margins(map, omega, mdot, Tt_in, pt_in)
    return (
        PR=root.PR,
        eta=root.eta,
        mdot=mdot,
        ht_out=ht_out,
        tau=tau,
        power=power,
        mdot_surge=margins.mdot_surge,
        mdot_choke=margins.mdot_choke,
        surge_margin=margins.surge_margin,
        choke_margin=margins.choke_margin,
        surge_margin_norm=margins.surge_margin_norm,
        choke_margin_norm=margins.choke_margin_norm,
        stall_flag=margins.stall_flag,
        choke_flag=margins.choke_flag,
    )
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
    n_points::Int=75,
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
        diagnostics = NamedTuple[]
        prior_roots = Float64[]
        for omega in omegas
            cond_diag = _compressor_condition_diagnostics(
                map,
                omega,
                Tt_in_f,
                pt_in_f,
                target_pr_f,
            )
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

            push!(
                diagnostics,
                (
                    omega=omega,
                    converged=(found.converged && !isempty(found.roots)),
                    failure_reason=String(found.failure_reason),
                    backoff_used=found.backoff_used,
                    requested_pr=found.target_pr,
                    requested_pt_out=found.target_pt_out,
                    achieved_pr=found.PR,
                    achieved_pt_out=found.pt_out,
                    n_target_roots=found.n_target_roots,
                    n_achieved_roots=length(found.roots),
                    flow_mdot_min=cond_diag.flow_mdot_min,
                    flow_mdot_max=cond_diag.flow_mdot_max,
                    pr_surge=cond_diag.pr_surge,
                    pr_choke=cond_diag.pr_choke,
                    pr_max=cond_diag.pr_max,
                    mdot_at_pr_max=cond_diag.mdot_at_pr_max,
                    pr_peak_boundary=cond_diag.pr_peak_boundary,
                    target_feasible_by_scan=cond_diag.target_feasible_by_scan,
                    root_bracketed_by_scan=cond_diag.root_bracketed_by_scan,
                    infeasibility_hint=cond_diag.infeasibility_hint,
                ),
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
                    mdot_surge=NaN,
                    mdot_choke=NaN,
                    surge_margin=NaN,
                    choke_margin=NaN,
                    surge_margin_norm=NaN,
                    choke_margin_norm=NaN,
                    stall_flag=false,
                    choke_flag=false,
                ))
                continue
            end

            selected = _compressor_select_roots(found.roots, :all)
            for (k, root) in enumerate(selected)
                out = _compressor_root_outputs(root, eos, map, pt_in_f, ht_in, omega, Tt_in_f)
                push!(rows, (
                    omega=omega,
                    branch_id=k,
                    PR=out.PR,
                    eta=out.eta,
                    mdot=out.mdot,
                    power=out.power,
                    converged=true,
                    backoff_used=found.backoff_used,
                    mdot_surge=out.mdot_surge,
                    mdot_choke=out.mdot_choke,
                    surge_margin=out.surge_margin,
                    choke_margin=out.choke_margin,
                    surge_margin_norm=out.surge_margin_norm,
                    choke_margin_norm=out.choke_margin_norm,
                    stall_flag=out.stall_flag,
                    choke_flag=out.choke_flag,
                ))
            end
            prior_roots = [r.mdot for r in found.roots]
        end
        return (mode=:all, branch=:all, rows=rows, diagnostics=diagnostics)
    end

    prs = fill(NaN, n_points)
    etas = fill(NaN, n_points)
    mdots = fill(NaN, n_points)
    powers = fill(NaN, n_points)
    converged = fill(false, n_points)
    backoff_used = fill(false, n_points)
    mdot_surge = fill(NaN, n_points)
    mdot_choke = fill(NaN, n_points)
    surge_margin = fill(NaN, n_points)
    choke_margin = fill(NaN, n_points)
    surge_margin_norm = fill(NaN, n_points)
    choke_margin_norm = fill(NaN, n_points)
    stall_flag = fill(false, n_points)
    choke_flag = fill(false, n_points)
    diagnostics = NamedTuple[]
    prior_roots = Float64[]

    for (i, omega) in enumerate(omegas)
        cond_diag = _compressor_condition_diagnostics(
            map,
            omega,
            Tt_in_f,
            pt_in_f,
            target_pr_f,
        )
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

        push!(
            diagnostics,
            (
                omega=omega,
                converged=(found.converged && !isempty(found.roots)),
                failure_reason=String(found.failure_reason),
                backoff_used=found.backoff_used,
                requested_pr=found.target_pr,
                requested_pt_out=found.target_pt_out,
                achieved_pr=found.PR,
                achieved_pt_out=found.pt_out,
                n_target_roots=found.n_target_roots,
                n_achieved_roots=length(found.roots),
                flow_mdot_min=cond_diag.flow_mdot_min,
                flow_mdot_max=cond_diag.flow_mdot_max,
                pr_surge=cond_diag.pr_surge,
                pr_choke=cond_diag.pr_choke,
                pr_max=cond_diag.pr_max,
                mdot_at_pr_max=cond_diag.mdot_at_pr_max,
                pr_peak_boundary=cond_diag.pr_peak_boundary,
                target_feasible_by_scan=cond_diag.target_feasible_by_scan,
                root_bracketed_by_scan=cond_diag.root_bracketed_by_scan,
                infeasibility_hint=cond_diag.infeasibility_hint,
            ),
        )

        if !found.converged || isempty(found.roots)
            continue
        end

        root = only(_compressor_select_roots(found.roots, branch))
        out = _compressor_root_outputs(root, eos, map, pt_in_f, ht_in, omega, Tt_in_f)
        prs[i] = out.PR
        etas[i] = out.eta
        mdots[i] = out.mdot
        powers[i] = out.power
        converged[i] = true
        backoff_used[i] = found.backoff_used
        mdot_surge[i] = out.mdot_surge
        mdot_choke[i] = out.mdot_choke
        surge_margin[i] = out.surge_margin
        choke_margin[i] = out.choke_margin
        surge_margin_norm[i] = out.surge_margin_norm
        choke_margin_norm[i] = out.choke_margin_norm
        stall_flag[i] = out.stall_flag
        choke_flag[i] = out.choke_flag
        prior_roots = [r.mdot for r in found.roots]
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
        mdot_surge=mdot_surge,
        mdot_choke=mdot_choke,
        surge_margin=surge_margin,
        choke_margin=choke_margin,
        surge_margin_norm=surge_margin_norm,
        choke_margin_norm=choke_margin_norm,
        stall_flag=stall_flag,
        choke_flag=choke_flag,
        diagnostics=diagnostics,
    )
end

function solve_compressor_operating_sweep(
    map::AbstractCompressorPerformanceMap;
    kwargs...,
)
    eos = Fluids.ideal_EOS()[:air]
    return solve_compressor_operating_sweep(map, eos; kwargs...)
end
