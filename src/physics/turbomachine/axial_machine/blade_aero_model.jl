"""
Unified blade aerodynamic closure.

Physical meaning of fields:
- `deviation_ref` [rad]:
  Baseline exit-flow deviation from the blade metal exit angle at near-design
  incidence. Positive values mean the flow exits at a smaller angle than the
  metal angle (`theta_out = theta_metal_out - deviation`).
- `deviation_incidence_sensitivity` [-]:
  Linear sensitivity of deviation to inlet incidence. With the current closure,
  `deviation = deviation_ref + deviation_incidence_sensitivity * incidence`.
- `loss_entropy_base` [J/(kg*K)]:
  Baseline stagnation-entropy rise across the row at near-design incidence.
- `loss_entropy_incidence` [J/(kg*K*rad^2)]:
  Quadratic off-design penalty on incidence in
  `delta_s_t = loss_entropy_base + loss_entropy_incidence * incidence^2`.
- `stall_incidence_limit` [rad]:
  Magnitude of incidence where stall margin reaches zero in this model.
  Reported stall margin is `stall_incidence_limit - abs(incidence)`.
- `theta_min` [rad], `theta_max` [rad]:
  Hard bounds on the modeled exit flow angle. These keep turning
  physically/numerically bounded near extreme incidence before conversion to
  `k_theta_exit = tan(theta_out)`.
"""
struct BladeAeroModel{T<:Real}
    deviation_ref::T
    deviation_incidence_sensitivity::T
    loss_entropy_base::T
    loss_entropy_incidence::T
    stall_incidence_limit::T
    theta_min::T
    theta_max::T
end

"""
Result of evaluating the blade aerodynamic closure at one row inlet state.

Field meanings:
- `theta_out` [rad]:
  Clamped exit flow angle used by the streamtube march.
- `theta_out_unclamped` [rad]:
  Unclamped exit flow angle predicted before applying the hard turning limits.
- `delta_s_t` [J/(kg*K)]:
  Modeled stagnation-entropy rise across the row.
- `stall_margin` [rad]:
  Signed incidence margin to the configured stall limit. Negative values mean
  the row is beyond the model's nominal stall incidence envelope.
- `regime`:
  Diagnostic flow regime selected by the closure. Currently `:normal` or
  `:stall`.
- `incidence` [rad]:
  Inlet relative flow angle minus metal inlet angle.
- `deviation` [rad]:
  Exit flow deviation from the metal exit angle.
- `theta_in` [rad]:
  Inlet relative flow angle seen by the blade row.
- `theta_metal_in`, `theta_metal_out` [rad]:
  Metal inlet/exit angles used to interpret incidence and deviation.
"""
struct BladeAeroResult{T<:Real}
    theta_out::T
    theta_out_unclamped::T
    delta_s_t::T
    stall_margin::T
    regime::Symbol
    incidence::T
    deviation::T
    theta_in::T
    theta_metal_in::T
    theta_metal_out::T
end

"""
Convenience constructor with rotor-like defaults.
"""
function rotor_aero_model(;
    deviation_ref::Real=0.0,
    deviation_incidence_sensitivity::Real=0.75,
    loss_entropy_base::Real=10.0,
    loss_entropy_incidence::Real=180.0,
    stall_incidence_limit::Real=0.32,
    theta_min::Real=deg2rad(-68.0),
    theta_max::Real=deg2rad(56.0),
)
    return BladeAeroModel{Float64}(
        Float64(deviation_ref),
        Float64(deviation_incidence_sensitivity),
        Float64(loss_entropy_base),
        Float64(loss_entropy_incidence),
        Float64(stall_incidence_limit),
        Float64(theta_min),
        Float64(theta_max),
    )
end

"""
Convenience constructor with stator-like defaults.
"""
function stator_aero_model(;
    deviation_ref::Real=0.0,
    deviation_incidence_sensitivity::Real=0.85,
    loss_entropy_base::Real=6.0,
    loss_entropy_incidence::Real=120.0,
    stall_incidence_limit::Real=0.30,
    theta_min::Real=deg2rad(-50.0),
    theta_max::Real=deg2rad(68.0),
)
    return BladeAeroModel{Float64}(
        Float64(deviation_ref),
        Float64(deviation_incidence_sensitivity),
        Float64(loss_entropy_base),
        Float64(loss_entropy_incidence),
        Float64(stall_incidence_limit),
        Float64(theta_min),
        Float64(theta_max),
    )
end

"""
Evaluate the row aerodynamic closure from dimensional inlet velocities and blade
metal angles.

This closure is intentionally simple:
- inlet flow angle is computed in the blade-relative frame
- incidence is measured against the metal inlet angle
- deviation is modeled as a linear function of incidence
- exit angle is clamped in angle space
- row loss is returned directly as a dimensional stagnation-entropy rise

Inputs:
- `Vx_in`, `Vtheta_in`: absolute inlet velocity components [m/s]
- `U`: local blade speed [m/s]

Returned quantities:
- `theta_out`: clamped exit flow angle used by the streamtube march
- `theta_out_unclamped`: raw deviation-model exit flow angle before clamping
- `delta_s_t`: stagnation-entropy rise across the row [J/(kg*K)]
- `stall_margin`: incidence margin relative to the configured stall limit
- `regime`: symbolic diagnostic regime (`:normal` or `:stall`)
- `incidence`, `deviation`, `theta_in`: angle diagnostics used by the closure
"""
function blade_aero(
    model::BladeAeroModel{T},
    theta_metal_in::U,
    theta_metal_out::U,
    Vx_in::U,
    Vtheta_in::U,
    U_blade::U,
) where {T<:Real,U<:Real}
    # Form the inlet relative-flow angle from the dimensional velocity triangle.
    theta_in = atan(Vtheta_in - U_blade, Vx_in)

    # Measure incidence from the metal inlet angle, then predict deviation from
    # the metal exit angle. Positive deviation reduces the flow exit angle
    # relative to the metal angle.
    incidence = theta_in - theta_metal_in
    stall_margin = model.stall_incidence_limit - abs(incidence)
    normal_deviation = model.deviation_ref + model.deviation_incidence_sensitivity * incidence
    normal_theta_out = theta_metal_out - normal_deviation
    normal_delta_s_t = model.loss_entropy_base + model.loss_entropy_incidence * incidence^2

    if stall_margin >= zero(stall_margin)
        regime = :normal
        theta_out_unclamped = normal_theta_out
        theta_out = clamp(theta_out_unclamped, model.theta_min, model.theta_max)
        delta_s_t = normal_delta_s_t
    else
        regime = :stall
        # In stalled flow, collapse turning toward the inlet flow direction and
        # increase entropy generation sharply as incidence exceeds the limit.
        overshoot = abs(incidence) - model.stall_incidence_limit
        overshoot_scale = max(abs(model.stall_incidence_limit), T(1e-6))
        severity = overshoot / overshoot_scale
        relax = inv(one(T) + T(3) * severity)
        theta_out_unclamped = theta_in + relax * (normal_theta_out - theta_in)
        theta_out = clamp(theta_out_unclamped, model.theta_min, model.theta_max)
        delta_s_t = max(
            normal_delta_s_t,
            (model.loss_entropy_base + model.loss_entropy_incidence * model.stall_incidence_limit^2) *
            (one(T) + T(8) * severity^2),
        )
    end
    deviation = theta_metal_out - theta_out
    return BladeAeroResult(
        theta_out,
        theta_out_unclamped,
        max(delta_s_t, zero(delta_s_t)),
        stall_margin,
        regime,
        incidence,
        deviation,
        theta_in,
        theta_metal_in,
        theta_metal_out,
    )
end
