"""
Axial-machine streamtube model definition.

Field meanings:
- `gamma`: ratio of specific heats for the working fluid (`cp/cv`), assumed
  constant across the solve.
- `gas_constant`: specific gas constant `R` [J/(kg*K)] for the working fluid.
- `A_ref`: reference annulus area [m^2] used to dimensionalize station areas.
- `A_station`: non-dimensional station area multipliers. Physical station area
  is `A_ref * A_station[k]` at station `k`. Length must be `n_rows + 1`.
- `rows`: ordered row models from inlet to outlet. Each row advances the
  solution from station `k` to `k+1`.
- `nu_theta_station1`: inlet circumferential velocity coefficient at station 1,
  with `nu_theta = V_theta / a0_in` (signed by swirl direction).
- `m_tip_bounds`: tabulation/operating domain bounds for reference tip Mach-like
  speed, `m_tip = omega * r_tip_ref / a0_in`.
- `phi_in_bounds`: tabulation/operating domain bounds for inlet flow
  coefficient, `phi_in = Vx_1 / U_m1`.
"""
struct AxialMachineModel
    # Thermodynamic properties for the quasi-1D flow model.
    gamma::Float64
    gas_constant::Float64

    # Geometry scaling: station area is A_ref * A_station[k].
    A_ref::Float64
    A_station::Vector{Float64}

    # Ordered aero-row closures that map station k -> k+1.
    rows::Vector{AxialRow}

    # Inlet swirl state (nu_theta) at station 1.
    nu_theta_station1::Float64

    # Recommended (m_tip, phi_in) bounds for map tabulation/sweeps.
    m_tip_bounds::Tuple{Float64,Float64}
    phi_in_bounds::Tuple{Float64,Float64}
end

function AxialMachineModel(
    gamma::Real,
    gas_constant::Real,
    A_ref::Real,
    A_station::Vector{<:Real},
    rows::Vector{AxialRow},
    nu_theta_station1::Real,
    m_tip_bounds::Tuple{<:Real,<:Real},
    phi_in_bounds::Tuple{<:Real,<:Real},
)
    gamma > 1 || error("gamma must be > 1")
    gas_constant > 0 || error("gas_constant must be > 0")
    A_ref > 0 || error("A_ref must be > 0")
    isempty(rows) && error("rows must not be empty")
    length(A_station) == length(rows) + 1 || error("A_station length must be n_rows + 1")
    all(a -> a > 0, A_station) || error("A_station entries must be > 0")
    m_lo, m_hi = Float64(m_tip_bounds[1]), Float64(m_tip_bounds[2])
    phi_lo, phi_hi = Float64(phi_in_bounds[1]), Float64(phi_in_bounds[2])
    m_hi > m_lo > 0 || error("m_tip_bounds must satisfy 0 < lo < hi")
    phi_hi > phi_lo > 0 || error("phi_in_bounds must satisfy 0 < lo < hi")
    return AxialMachineModel(
        Float64(gamma),
        Float64(gas_constant),
        Float64(A_ref),
        Float64.(A_station),
        rows,
        Float64(nu_theta_station1),
        (m_lo, m_hi),
        (phi_lo, phi_hi),
    )
end

function first_rotor_index(model::AxialMachineModel)
    idx = findfirst(row -> row.kind == :rotor, model.rows)
    return isnothing(idx) ? 1 : idx
end
