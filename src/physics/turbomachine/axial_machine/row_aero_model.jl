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
Row aerodynamic input state.
"""
struct RowAeroInput{T}
    nu_x_in::T
    nu_theta_in::T
    nu_u::T
    tau_in::T
    pi_in::T
    mu::T
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

"""
Axial row descriptor.
"""
struct AxialRow{M<:AbstractRowAeroModel}
    kind::Symbol
    aero::M
    r_mean::Float64
    r_tip::Float64
    omega_sign::Int8
end

function AxialRow(
    kind::Symbol,
    aero::AbstractRowAeroModel,
    r_mean::Real,
    r_tip::Real,
    omega_sign::Integer,
)
    kind in (:rotor, :stator) || error("row kind must be :rotor or :stator")
    r_mean > 0 || error("row r_mean must be > 0")
    r_tip > 0 || error("row r_tip must be > 0")
    sign_i8 = Int8(omega_sign)
    if kind == :rotor
        sign_i8 in Int8[-1, 1] || error("rotor omega_sign must be -1 or +1")
        aero isa RotorAeroModel || error("rotor rows require RotorAeroModel")
    else
        sign_i8 == 0 || error("stator omega_sign must be 0")
        aero isa StatorAeroModel || error("stator rows require StatorAeroModel")
    end
    return AxialRow{typeof(aero)}(kind, aero, Float64(r_mean), Float64(r_tip), sign_i8)
end

function row_aero(
    model::RotorAeroModel{T},
    row::AxialRow,
    input::RowAeroInput{U},
) where {T<:Real,U<:Real}
    beta_in = atan(input.nu_theta_in - input.nu_u, input.nu_x_in)
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
    row::AxialRow,
    input::RowAeroInput{U},
) where {T<:Real,U<:Real}
    alpha_in = atan(input.nu_theta_in, input.nu_x_in)
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
