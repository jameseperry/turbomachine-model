"""
Axial row descriptor.

Field meanings:
- `kind`: row category (`:rotor` or `:stator`) used for validation and
  interpretation of kinematics.
- `aero`: aerodynamic closure model for this row (`BladeAeroModel`) that maps
  inlet flow state to turning/loss outputs.
- `r_hub`: hub radius [m] of the annulus represented by this row.
- `r_tip`: tip radius [m] for this row.
- `speed_ratio_to_ref`: row shaft speed divided by the reference shaft speed
  (`omega_row / omega_ref`). Use `0.0` for stators, negative values for
  counter-rotation, and values with magnitude > 1 for geared-up rows.
"""
struct AxialRow
    kind::Symbol
    aero::BladeAeroModel
    r_hub::Float64
    r_tip::Float64
    speed_ratio_to_ref::Float64
end

function AxialRow(
    kind::Symbol,
    aero::BladeAeroModel,
    r_hub::Real,
    r_tip::Real,
    speed_ratio_to_ref::Real,
)
    kind in (:rotor, :stator) || error("row kind must be :rotor or :stator")
    r_hub >= 0 || error("row r_hub must be >= 0")
    r_tip > 0 || error("row r_tip must be > 0")
    r_tip > r_hub || error("row r_tip must be > r_hub")
    speed_ratio = Float64(speed_ratio_to_ref)
    if kind == :rotor
        speed_ratio != 0.0 || error("rotor speed_ratio_to_ref must be nonzero")
    else
        speed_ratio == 0.0 || error("stator speed_ratio_to_ref must be 0.0")
    end
    return AxialRow(
        kind,
        aero,
        Float64(r_hub),
        Float64(r_tip),
        speed_ratio,
    )
end

"""
Axial-machine streamtube model definition.

Field meanings:
- `gamma`: ratio of specific heats for the working fluid (`cp/cv`), assumed
  constant across the solve.
- `gas_constant`: specific gas constant `R` [J/(kg*K)] for the working fluid.
- `r_tip_ref`: reference rotor-tip radius [m] used for non-dimensional speed
  scaling (`m_tip = omega_ref * r_tip_ref / a0_in`).
- `rows`: ordered row models from inlet to outlet. Each row advances the
  solution from station `k` to `k+1`.
- `first_rotor_index`: cached index of the first rotor row. If no rotor row is
  present, this falls back to `1` for compatibility with existing behavior.
- `m_tip_bounds`: tabulation/operating domain bounds for reference tip Mach-like
  speed.
- `phi_in_bounds`: tabulation/operating domain bounds for inlet flow
  coefficient.
"""
struct AxialMachineModel
    gamma::Float64
    gas_constant::Float64
    r_tip_ref::Float64
    rows::Vector{AxialRow}
    first_rotor_index::Int
    m_tip_bounds::Tuple{Float64,Float64}
    phi_in_bounds::Tuple{Float64,Float64}
end

function AxialMachineModel(
    gamma::Real,
    gas_constant::Real,
    r_tip_ref::Real,
    rows::Vector{AxialRow},
    m_tip_bounds::Tuple{<:Real,<:Real},
    phi_in_bounds::Tuple{<:Real,<:Real},
)
    gamma > 1 || error("gamma must be > 1")
    gas_constant > 0 || error("gas_constant must be > 0")
    r_tip_ref > 0 || error("r_tip_ref must be > 0")
    isempty(rows) && error("rows must not be empty")
    m_lo, m_hi = Float64(m_tip_bounds[1]), Float64(m_tip_bounds[2])
    phi_lo, phi_hi = Float64(phi_in_bounds[1]), Float64(phi_in_bounds[2])
    m_hi > m_lo > 0 || error("m_tip_bounds must satisfy 0 < lo < hi")
    phi_hi > phi_lo > 0 || error("phi_in_bounds must satisfy 0 < lo < hi")
    first_rotor_idx = let idx = findfirst(row -> row.kind == :rotor, rows)
        isnothing(idx) ? 1 : idx
    end
    return AxialMachineModel(
        Float64(gamma),
        Float64(gas_constant),
        Float64(r_tip_ref),
        rows,
        first_rotor_idx,
        (m_lo, m_hi),
        (phi_lo, phi_hi),
    )
end

function meanline_radii(model::AxialMachineModel)
    return [0.5 * (row.r_hub + row.r_tip) for row in model.rows]
end

row_annulus_area(row::AxialRow) = pi * (row.r_tip^2 - row.r_hub^2)

function station_area(model::AxialMachineModel, station_index::Integer)
    n_rows = length(model.rows)
    1 <= station_index <= (n_rows + 1) || error("station_index out of bounds")
    if station_index == 1
        return row_annulus_area(model.rows[1])
    elseif station_index == n_rows + 1
        return row_annulus_area(model.rows[end])
    end
    a_prev = row_annulus_area(model.rows[station_index - 1])
    a_next = row_annulus_area(model.rows[station_index])
    return 0.5 * (a_prev + a_next)
end
