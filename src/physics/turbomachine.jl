"""
Turbomachine residual equations.
"""

"""
Compute turbomachine residuals.

Returns `(R_pr, R_eta, R_P)`:
- `R_pr  = pt_out - PR * pt_in`
- `R_eta = ht_out - ht_in - (h2s - ht_in) / eta`
- `R_P   = tau * omega - mdot * (ht_out - ht_in)`
"""
function turbomachine_residuals(
    performance_map::Fluids.PerformanceMap,
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
    map_vals = Fluids.map_pr_eta_from_stagnation(
        performance_map,
        omega,
        mdot,
        Tt_in,
        pt_in,
    )
    PR = map_vals.PR
    eta = map_vals.eta
    eta > 0 || error("map eta must be > 0")

    h2s = Fluids.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)

    R_pr = pt_out - PR * pt_in
    R_eta = ht_out - ht_in - (h2s - ht_in) / eta
    R_P = tau * omega - mdot * (ht_out - ht_in)
    return (R_pr, R_eta, R_P)
end
