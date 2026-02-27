"""
Compressor residual equations.
"""

@inline _primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x

function compressor_residual_scales(
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
Compute compressor residuals.

Returns `(R_pout, R_dh_eff, R_P)`:
- `R_pout = pt_out - PR * pt_in`
- `R_dh_eff = eta * (ht_out - ht_in) - (h2s - ht_in)`
- `R_P = tau * omega - mdot * (ht_out - ht_in)`
"""
function compressor_residuals(
    compressor_map::AbstractCompressorPerformanceMap,
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
    map_vals = compressor_performance_map_from_stagnation(
        compressor_map,
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
Compute scaled compressor residuals.
"""
function compressor_residuals_scaled(
    compressor_map::AbstractCompressorPerformanceMap,
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
    scales = compressor_residual_scales(
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
