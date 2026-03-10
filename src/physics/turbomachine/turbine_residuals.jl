"""
Turbine residual equations.
"""

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
    turbine_map::TabulatedPerformanceMap,
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
    map_vals = performance_from_stagnation(
        turbine_map,
        omega,
        mdot,
        Tt_in,
        pt_in,
    )
    eta = map_vals.eta
    _primal_value(eta) > 0 || error("map eta must be > 0")

    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)

    R_pout = pt_out - map_vals.pressure_ratio * pt_in
    R_dh_eff = (ht_in - ht_out) - eta * (ht_in - h2s)
    R_P = tau * omega - mdot * (ht_out - ht_in)
    return (R_pout, R_dh_eff, R_P)
end

"""
Compute scaled turbine residuals.
"""
function turbine_residuals_scaled(
    turbine_map::TabulatedPerformanceMap,
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
    R_pout, R_dh_eff, R_P = turbine_residuals(
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
        R_pout / scales.massflow_scale,
        R_dh_eff / scales.enthalpy_scale,
        R_P / scales.power_scale,
    )
end
