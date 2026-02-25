"""
Turbomachine residual equations.
"""

using NonlinearSolve

@inline _primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x

function turbomachine_residual_scales(
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    ht_out::Real,
    mdot::Real,
    omega::Real,
    tau::Real;
    pressure_scale::Union{Nothing,Real}=nothing,
    enthalpy_scale::Union{Nothing,Real}=nothing,
    power_scale::Union{Nothing,Real}=nothing,
)
    p_ref = if pressure_scale === nothing
        max(abs(_primal_value(pt_in)), abs(_primal_value(pt_out)), 1.0)
    else
        _primal_value(pressure_scale) > 0 || error("pressure_scale must be > 0")
        pressure_scale
    end

    h_ref = if enthalpy_scale === nothing
        max(abs(_primal_value(ht_in)), abs(_primal_value(ht_out)), 1.0)
    else
        _primal_value(enthalpy_scale) > 0 || error("enthalpy_scale must be > 0")
        enthalpy_scale
    end

    power_in = _primal_value(tau * omega)
    power_out = _primal_value(mdot * (ht_out - ht_in))
    P_ref = if power_scale === nothing
        max(abs(power_in), abs(power_out), 1.0)
    else
        _primal_value(power_scale) > 0 || error("power_scale must be > 0")
        power_scale
    end

    return (pressure_scale=p_ref, enthalpy_scale=h_ref, power_scale=P_ref)
end

"""
Compute turbomachine residuals.

Returns `(R_pout, R_dh_eff, R_P)`:
- `R_pout = pt_out - PR * pt_in`
- `R_dh_eff = eta * (ht_out - ht_in) - (h2s - ht_in)`
- `R_P   = tau * omega - mdot * (ht_out - ht_in)`
"""
function turbomachine_residuals(
    performance_map::AbstractPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    ht_out::Real,
    mdot::Real,
    omega::Real,
    tau::Real,
)
    Tt_in = Fluids.temperature(eos, pt_in, ht_in)
    map_vals = performance_map_from_stagnation(
        performance_map,
        omega,
        mdot,
        Tt_in,
        pt_in,
    )
    PR = map_vals.PR
    eta = map_vals.eta
    _primal_value(eta) > 0 || error("map eta must be > 0")

    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)

    R_pout = pt_out - PR * pt_in
    R_dh_eff = eta * (ht_out - ht_in) - (h2s - ht_in)
    R_P = tau * omega - mdot * (ht_out - ht_in)
    return (R_pout, R_dh_eff, R_P)
end

"""
Compute scaled turbomachine residuals.

Returns `(r_pout, r_dh_eff, r_P)` where:
- `r_pout = R_pout / p_ref`
- `r_dh_eff = R_dh_eff / h_ref`
- `r_P    = R_P / P_ref`

Default scale factors are derived from the local state:
- `p_ref = max(abs(pt_in), abs(pt_out), 1.0)`
- `h_ref = max(abs(ht_in), abs(ht_out), 1.0)`
- `P_ref = max(abs(tau * omega), abs(mdot * (ht_out - ht_in)), 1.0)`

You can override each scale with keyword arguments:
- `pressure_scale`
- `enthalpy_scale`
- `power_scale`
"""
function turbomachine_residuals_scaled(
    performance_map::AbstractPerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    ht_out::Real,
    mdot::Real,
    omega::Real,
    tau::Real;
    pressure_scale::Union{Nothing,Real}=nothing,
    enthalpy_scale::Union{Nothing,Real}=nothing,
    power_scale::Union{Nothing,Real}=nothing,
)
    R_pout, R_dh_eff, R_P = turbomachine_residuals(
        performance_map,
        eos,
        pt_in,
        ht_in,
        pt_out,
        ht_out,
        mdot,
        omega,
        tau,
    )
    scales = turbomachine_residual_scales(
        pt_in,
        ht_in,
        pt_out,
        ht_out,
        mdot,
        omega,
        tau;
        pressure_scale=pressure_scale,
        enthalpy_scale=enthalpy_scale,
        power_scale=power_scale,
    )
    return (
        R_pout / scales.pressure_scale,
        R_dh_eff / scales.enthalpy_scale,
        R_P / scales.power_scale,
    )
end

"""
Solve a simple turbomachine operating point with four boundary conditions.

Fixed boundary conditions:
- `pt_in`
- `ht_in`
- `pt_out`
- `omega`

Solved unknowns:
- `mdot`
- `ht_out`
- `tau`

By default this solve uses scaled residuals (`scaled_residuals=true`).
You can disable scaling or provide explicit residual scales through:
- `scaled_residuals`
- `pressure_scale`
- `enthalpy_scale`
- `power_scale`
"""
function solve_turbomachine_operating_point(
    performance_map::AbstractPerformanceMap,
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
    scale_overrides = turbomachine_residual_scales(
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
            r_pout, r_dh_eff, r_P = turbomachine_residuals_scaled(
                performance_map,
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
            R_pout, R_dh_eff, R_P = turbomachine_residuals(
                performance_map,
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

    residuals = turbomachine_residuals(
        performance_map,
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
    map_vals = performance_map_from_stagnation(
        performance_map,
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
