"""
Abstract aerodynamic row model used by axial streamtube marching.
"""
abstract type AbstractRowAeroModel end

"""
Simple rotor aerodynamic closure.
"""
Base.@kwdef struct RotorAeroModel{T<:Real} <: AbstractRowAeroModel
    beta_ref::T = T(-0.55)
    beta_incidence_sensitivity::T = T(0.75)
    loss_base::T = T(0.010)
    loss_incidence::T = T(0.18)
    stall_incidence_limit::T = T(0.32)
    k_theta_min::T = T(-2.5)
    k_theta_max::T = T(1.5)
end

"""
Simple stator aerodynamic closure.
"""
Base.@kwdef struct StatorAeroModel{T<:Real} <: AbstractRowAeroModel
    alpha_ref::T = T(0.45)
    alpha_incidence_sensitivity::T = T(0.85)
    loss_base::T = T(0.006)
    loss_incidence::T = T(0.12)
    stall_incidence_limit::T = T(0.30)
    k_theta_min::T = T(-1.2)
    k_theta_max::T = T(2.5)
end

"""
Row aerodynamic output contract for the marcher.
"""
struct RowAeroOutput{T}
    k_theta_exit::T
    delta_s_hat::T
    stall_margin::T
    valid::Bool
    diagnostics::NamedTuple
end

function row_aero(
    model::RotorAeroModel{T},
    nu_x_in::U,
    nu_theta_in::U,
    nu_u::U,
) where {T<:Real,U<:Real}
    beta_in = atan(nu_theta_in - nu_u, nu_x_in)
    incidence = beta_in - model.beta_ref
    beta_out = model.beta_ref - model.beta_incidence_sensitivity * incidence
    k_theta = clamp(tan(beta_out), model.k_theta_min, model.k_theta_max)
    delta_s_hat = model.loss_base + model.loss_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(k_theta) && isfinite(delta_s_hat)
    return RowAeroOutput(
        k_theta,
        max(delta_s_hat, zero(delta_s_hat)),
        stall_margin,
        valid,
        (incidence=incidence, beta_in=beta_in, beta_out=beta_out),
    )
end

function row_aero(
    model::StatorAeroModel{T},
    nu_x_in::U,
    nu_theta_in::U,
    _nu_u::U,
) where {T<:Real,U<:Real}
    alpha_in = atan(nu_theta_in, nu_x_in)
    incidence = alpha_in - model.alpha_ref
    alpha_out = model.alpha_ref - model.alpha_incidence_sensitivity * incidence
    k_theta = clamp(tan(alpha_out), model.k_theta_min, model.k_theta_max)
    delta_s_hat = model.loss_base + model.loss_incidence * incidence^2
    stall_margin = model.stall_incidence_limit - abs(incidence)
    valid = isfinite(k_theta) && isfinite(delta_s_hat)
    return RowAeroOutput(
        k_theta,
        max(delta_s_hat, zero(delta_s_hat)),
        stall_margin,
        valid,
        (incidence=incidence, alpha_in=alpha_in, alpha_out=alpha_out),
    )
end
