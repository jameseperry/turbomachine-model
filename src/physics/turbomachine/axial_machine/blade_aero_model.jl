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
- `loss_base` [-]:
  Baseline non-dimensional entropy-loss term at near-design incidence.
- `loss_incidence` [-]:
  Quadratic penalty coefficient on off-design incidence in
  `delta_s_hat = loss_base + loss_incidence * incidence^2`.
- `stall_incidence_limit` [rad]:
  Magnitude of incidence where stall margin reaches zero in this model.
  Reported stall margin is `stall_incidence_limit - abs(incidence)`.
- `k_theta_min` [-], `k_theta_max` [-]:
  Hard bounds on `k_theta_exit = tan(theta_out)` used to keep turning
  physically/numerically bounded near extreme incidence.
"""
struct BladeAeroModel{T<:Real}
    deviation_ref::T
    deviation_incidence_sensitivity::T
    loss_base::T
    loss_incidence::T
    stall_incidence_limit::T
    k_theta_min::T
    k_theta_max::T
end

"""
Convenience constructor with rotor-like defaults.
"""
function rotor_aero_model(;
    deviation_ref::Real=0.0,
    deviation_incidence_sensitivity::Real=0.75,
    loss_base::Real=0.010,
    loss_incidence::Real=0.18,
    stall_incidence_limit::Real=0.32,
    k_theta_min::Real=-2.5,
    k_theta_max::Real=1.5,
)
    return BladeAeroModel{Float64}(
        Float64(deviation_ref),
        Float64(deviation_incidence_sensitivity),
        Float64(loss_base),
        Float64(loss_incidence),
        Float64(stall_incidence_limit),
        Float64(k_theta_min),
        Float64(k_theta_max),
    )
end

"""
Convenience constructor with stator-like defaults.
"""
function stator_aero_model(;
    deviation_ref::Real=0.0,
    deviation_incidence_sensitivity::Real=0.85,
    loss_base::Real=0.006,
    loss_incidence::Real=0.12,
    stall_incidence_limit::Real=0.30,
    k_theta_min::Real=-1.2,
    k_theta_max::Real=2.5,
)
    return BladeAeroModel{Float64}(
        Float64(deviation_ref),
        Float64(deviation_incidence_sensitivity),
        Float64(loss_base),
        Float64(loss_incidence),
        Float64(stall_incidence_limit),
        Float64(k_theta_min),
        Float64(k_theta_max),
    )
end

"""
Evaluate blade aerodynamics from incoming velocity components and blade metal
angles.
"""
function blade_aero(
    model::BladeAeroModel{T},
    theta_metal_in::U,
    theta_metal_out::U,
    nu_x_in::U,
    nu_theta_in::U,
    nu_u::U,
) where {T<:Real,U<:Real}
    theta_in = atan(nu_theta_in - nu_u, nu_x_in)
    incidence = theta_in - theta_metal_in
    deviation = model.deviation_ref + model.deviation_incidence_sensitivity * incidence
    theta_out = theta_metal_out - deviation
    k_theta = clamp(tan(theta_out), model.k_theta_min, model.k_theta_max)
    delta_s_hat = model.loss_base + model.loss_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(k_theta) && isfinite(delta_s_hat)
    return (
        k_theta_exit=k_theta,
        delta_s_hat=max(delta_s_hat, zero(delta_s_hat)),
        stall_margin=stall_margin,
        valid=valid,
        diagnostics=(
            incidence=incidence,
            deviation=deviation,
            theta_in=theta_in,
            theta_out=theta_out,
            theta_metal_in=theta_metal_in,
            theta_metal_out=theta_metal_out,
        ),
    )
end
