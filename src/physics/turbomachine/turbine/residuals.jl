"""
Turbine residual equations and operating-point solve helpers.
"""

using NonlinearSolve

@inline _primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x

function turbine_residual_scales(
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    ht_out::Real,
    mdot::Real,
    omega::Real,
    tau::Real;
    massflow_scale::Union{Nothing,Real}=nothing,
    enthalpy_scale::Union{Nothing,Real}=nothing,
    power_scale::Union{Nothing,Real}=nothing,
)
    mdot_ref = if massflow_scale === nothing
        max(abs(_primal_value(mdot)), 1.0)
    else
        _primal_value(massflow_scale) > 0 || error("massflow_scale must be > 0")
        massflow_scale
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

    return (massflow_scale=mdot_ref, enthalpy_scale=h_ref, power_scale=P_ref)
end

"""
Compute turbine residuals.

Returns `(R_mdot_map, R_dh_eff, R_P)`:
- `R_mdot_map = mdot - mdot_map`
- `R_dh_eff = (ht_in - ht_out) - eta * (ht_in - h2s)`
- `R_P = tau * omega - mdot * (ht_out - ht_in)`
"""
function turbine_residuals(
    turbine_map::AbstractTurbinePerformanceMap,
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
    map_vals = turbine_performance_map_from_stagnation(
        turbine_map,
        omega,
        pt_in,
        pt_out,
        Tt_in,
    )
    eta = map_vals.eta
    _primal_value(eta) > 0 || error("map eta must be > 0")

    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)

    R_mdot_map = mdot - map_vals.mdot
    R_dh_eff = (ht_in - ht_out) - eta * (ht_in - h2s)
    R_P = tau * omega - mdot * (ht_out - ht_in)
    return (R_mdot_map, R_dh_eff, R_P)
end

"""
Compute scaled turbine residuals.
"""
function turbine_residuals_scaled(
    turbine_map::AbstractTurbinePerformanceMap,
    eos::Fluids.AbstractEOS,
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    ht_out::Real,
    mdot::Real,
    omega::Real,
    tau::Real;
    massflow_scale::Union{Nothing,Real}=nothing,
    enthalpy_scale::Union{Nothing,Real}=nothing,
    power_scale::Union{Nothing,Real}=nothing,
)
    R_mdot_map, R_dh_eff, R_P = turbine_residuals(
        turbine_map,
        eos,
        pt_in,
        ht_in,
        pt_out,
        ht_out,
        mdot,
        omega,
        tau,
    )
    scales = turbine_residual_scales(
        pt_in,
        ht_in,
        pt_out,
        ht_out,
        mdot,
        omega,
        tau;
        massflow_scale=massflow_scale,
        enthalpy_scale=enthalpy_scale,
        power_scale=power_scale,
    )
    return (
        R_mdot_map / scales.massflow_scale,
        R_dh_eff / scales.enthalpy_scale,
        R_P / scales.power_scale,
    )
end

"""
Solve a simple turbine operating point with four boundary conditions.

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
function solve_turbine_operating_point(
    turbine_map::AbstractTurbinePerformanceMap,
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
    massflow_scale::Union{Nothing,Real}=nothing,
    enthalpy_scale::Union{Nothing,Real}=nothing,
    power_scale::Union{Nothing,Real}=nothing,
)
    scale_overrides = turbine_residual_scales(
        pt_in,
        ht_in,
        pt_out,
        ht_out_guess,
        mdot_guess,
        omega,
        tau_guess;
        massflow_scale=massflow_scale,
        enthalpy_scale=enthalpy_scale,
        power_scale=power_scale,
    )

    function f!(R, u, _)
        mdot = u[1]
        ht_out = u[2]
        tau = u[3]
        if scaled_residuals
            r_mdot_map, r_dh_eff, r_P = turbine_residuals_scaled(
                turbine_map,
                eos,
                pt_in,
                ht_in,
                pt_out,
                ht_out,
                mdot,
                omega,
                tau;
                massflow_scale=scale_overrides.massflow_scale,
                enthalpy_scale=scale_overrides.enthalpy_scale,
                power_scale=scale_overrides.power_scale,
            )
            R[1] = r_mdot_map
            R[2] = r_dh_eff
            R[3] = r_P
        else
            R_mdot_map, R_dh_eff, R_P = turbine_residuals(
                turbine_map,
                eos,
                pt_in,
                ht_in,
                pt_out,
                ht_out,
                mdot,
                omega,
                tau,
            )
            R[1] = R_mdot_map
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

    residuals = turbine_residuals(
        turbine_map,
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
    map_vals = turbine_performance_map_from_stagnation(
        turbine_map,
        omega,
        pt_in,
        pt_out,
        Tt_in,
    )

    return (
        mdot=mdot,
        ht_out=ht_out,
        tau=tau,
        PR_turb=map_vals.PR_turb,
        eta=map_vals.eta,
        residuals=residuals,
        retcode=sol.retcode,
        converged=(string(sol.retcode) == "Success"),
        solution=sol,
    )
end
