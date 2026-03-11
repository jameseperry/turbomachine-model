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
- `delta_s_t`: stagnation-entropy rise across the row [J/(kg*K)]
- `stall_margin`: incidence margin relative to the configured stall limit
- `diagnostics`: angles and incidence/deviation terms used by the closure
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
    deviation = model.deviation_ref + model.deviation_incidence_sensitivity * incidence
    theta_out_unclamped = theta_metal_out - deviation

    # Clamp the exit flow angle directly in angle space.
    theta_out = clamp(theta_out_unclamped, model.theta_min, model.theta_max)

    # Return a dimensional stagnation-entropy rise directly. The streamtube
    # march can add this to the row inlet stagnation entropy without any cp-based
    # nondimensional bridge.
    delta_s_t = model.loss_entropy_base + model.loss_entropy_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(theta_out) && isfinite(delta_s_t)
    return (
        theta_out=theta_out,
        delta_s_t=max(delta_s_t, zero(delta_s_t)),
        stall_margin=stall_margin,
        valid=valid,
        diagnostics=(
            incidence=incidence,
            deviation=deviation,
            theta_in=theta_in,
            theta_out=theta_out,
            theta_out_unclamped=theta_out_unclamped,
            theta_metal_in=theta_metal_in,
            theta_metal_out=theta_metal_out,
        ),
    )
end
