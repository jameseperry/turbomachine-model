"""
Unified blade aerodynamic closure.

Physical meaning of fields:
- `theta_ref` [rad]:
  Design relative inlet flow angle seen by the blade row. Incidence is
  computed as `theta_in - theta_ref`.
- `theta_incidence_sensitivity` [-]:
  Linear turning response to incidence. With the current closure,
  `theta_out = theta_ref - theta_incidence_sensitivity * incidence`.
  Larger values make exit angle (and therefore loading/turning) change more
  aggressively as incidence moves off design.
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
    theta_ref::T
    theta_incidence_sensitivity::T
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
    theta_ref::Real=-0.55,
    theta_incidence_sensitivity::Real=0.75,
    loss_base::Real=0.010,
    loss_incidence::Real=0.18,
    stall_incidence_limit::Real=0.32,
    k_theta_min::Real=-2.5,
    k_theta_max::Real=1.5,
)
    return BladeAeroModel{Float64}(
        Float64(theta_ref),
        Float64(theta_incidence_sensitivity),
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
    theta_ref::Real=0.45,
    theta_incidence_sensitivity::Real=0.85,
    loss_base::Real=0.006,
    loss_incidence::Real=0.12,
    stall_incidence_limit::Real=0.30,
    k_theta_min::Real=-1.2,
    k_theta_max::Real=2.5,
)
    return BladeAeroModel{Float64}(
        Float64(theta_ref),
        Float64(theta_incidence_sensitivity),
        Float64(loss_base),
        Float64(loss_incidence),
        Float64(stall_incidence_limit),
        Float64(k_theta_min),
        Float64(k_theta_max),
    )
end

"""
Evaluate blade aerodynamics from incoming velocity components.
"""
function blade_aero(
    model::BladeAeroModel{T},
    nu_x_in::U,
    nu_theta_in::U,
    nu_u::U,
) where {T<:Real,U<:Real}
    theta_in = atan(nu_theta_in - nu_u, nu_x_in)
    incidence = theta_in - model.theta_ref
    theta_out = model.theta_ref - model.theta_incidence_sensitivity * incidence
    k_theta = clamp(tan(theta_out), model.k_theta_min, model.k_theta_max)
    delta_s_hat = model.loss_base + model.loss_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(k_theta) && isfinite(delta_s_hat)
    return (
        k_theta_exit=k_theta,
        delta_s_hat=max(delta_s_hat, zero(delta_s_hat)),
        stall_margin=stall_margin,
        valid=valid,
        diagnostics=(incidence=incidence, theta_in=theta_in, theta_out=theta_out),
    )
end
