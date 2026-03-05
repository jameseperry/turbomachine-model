using TOML
import ....Utility: write_toml, read_toml

function _aero_to_toml_dict(model::RotorAeroModel)
    return Dict{String,Any}(
        "format" => "rotor_aero_model",
        "beta_ref" => Float64(model.beta_ref),
        "beta_incidence_sensitivity" => Float64(model.beta_incidence_sensitivity),
        "loss_base" => Float64(model.loss_base),
        "loss_incidence" => Float64(model.loss_incidence),
        "stall_incidence_limit" => Float64(model.stall_incidence_limit),
        "k_theta_min" => Float64(model.k_theta_min),
        "k_theta_max" => Float64(model.k_theta_max),
    )
end

function _aero_to_toml_dict(model::StatorAeroModel)
    return Dict{String,Any}(
        "format" => "stator_aero_model",
        "alpha_ref" => Float64(model.alpha_ref),
        "alpha_incidence_sensitivity" => Float64(model.alpha_incidence_sensitivity),
        "loss_base" => Float64(model.loss_base),
        "loss_incidence" => Float64(model.loss_incidence),
        "stall_incidence_limit" => Float64(model.stall_incidence_limit),
        "k_theta_min" => Float64(model.k_theta_min),
        "k_theta_max" => Float64(model.k_theta_max),
    )
end

function _aero_from_toml_dict(data::Dict{String,Any})
    haskey(data, "format") || error("aero model missing format")
    fmt = String(data["format"])
    if fmt == "rotor_aero_model"
        return RotorAeroModel{Float64}(
            beta_ref=Float64(data["beta_ref"]),
            beta_incidence_sensitivity=Float64(data["beta_incidence_sensitivity"]),
            loss_base=Float64(data["loss_base"]),
            loss_incidence=Float64(data["loss_incidence"]),
            stall_incidence_limit=Float64(data["stall_incidence_limit"]),
            k_theta_min=Float64(data["k_theta_min"]),
            k_theta_max=Float64(data["k_theta_max"]),
        )
    elseif fmt == "stator_aero_model"
        return StatorAeroModel{Float64}(
            alpha_ref=Float64(data["alpha_ref"]),
            alpha_incidence_sensitivity=Float64(data["alpha_incidence_sensitivity"]),
            loss_base=Float64(data["loss_base"]),
            loss_incidence=Float64(data["loss_incidence"]),
            stall_incidence_limit=Float64(data["stall_incidence_limit"]),
            k_theta_min=Float64(data["k_theta_min"]),
            k_theta_max=Float64(data["k_theta_max"]),
        )
    end
    error("unsupported aero model format=$(fmt)")
end

function _row_to_toml_dict(row::AxialRow)
    return Dict{String,Any}(
        "kind" => String(row.kind),
        "r_mean" => row.r_mean,
        "r_tip" => row.r_tip,
        "omega_sign" => Int(row.omega_sign),
        "aero" => _aero_to_toml_dict(row.aero),
    )
end

function _row_from_toml_dict(data::Dict{String,Any})
    haskey(data, "kind") || error("row missing kind")
    haskey(data, "r_mean") || error("row missing r_mean")
    haskey(data, "r_tip") || error("row missing r_tip")
    haskey(data, "omega_sign") || error("row missing omega_sign")
    haskey(data, "aero") || error("row missing aero")
    kind = Symbol(String(data["kind"]))
    aero = _aero_from_toml_dict(data["aero"])
    return AxialRow(
        kind,
        aero,
        Float64(data["r_mean"]),
        Float64(data["r_tip"]),
        Int(data["omega_sign"]),
    )
end

function _find_or_create_axial_group!(data::Dict{String,Any}, group::AbstractString)
    isempty(group) && return data
    node = data
    for key in split(group, '.')
        if !haskey(node, key)
            node[key] = Dict{String,Any}()
        end
        child = node[key]
        child isa Dict || error("group path conflicts with non-table key $(key)")
        node = child
    end
    return node
end

function _find_axial_group(data::Dict{String,Any}, group::AbstractString)
    isempty(group) && return data
    node = data
    for key in split(group, '.')
        haskey(node, key) || error("missing TOML group $(group)")
        child = node[key]
        child isa Dict || error("TOML group $(group) is not a table")
        node = child
    end
    return node
end

function write_toml(
    model::AxialMachineModel,
    path::AbstractString;
    group::AbstractString="compressor_meanline_model",
)
    data = Dict{String,Any}()
    node = _find_or_create_axial_group!(data, group)
    # Keep format stable for backward compatibility with existing files/scripts.
    node["format"] = "compressor_meanline_model"
    node["format_version"] = 1
    node["gamma"] = model.gamma
    node["gas_constant"] = model.gas_constant
    node["A_ref"] = model.A_ref
    node["A_station"] = model.A_station
    node["nu_theta_first_rotor"] = model.nu_theta_first_rotor
    node["m_tip_bounds"] = [model.m_tip_bounds[1], model.m_tip_bounds[2]]
    node["phi_in_bounds"] = [model.phi_in_bounds[1], model.phi_in_bounds[2]]
    node["rows"] = [_row_to_toml_dict(row) for row in model.rows]
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{AxialMachineModel},
    path::AbstractString;
    group::AbstractString="compressor_meanline_model",
)
    data = TOML.parsefile(path)
    node = _find_axial_group(data, group)
    for key in ("gamma", "gas_constant", "A_ref", "A_station", "rows")
        haskey(node, key) || error("missing TOML key $(key)")
    end
    nu_theta_first_rotor = if haskey(node, "nu_theta_first_rotor")
        Float64(node["nu_theta_first_rotor"])
    elseif haskey(node, "nu_theta_station1")
        # Backward-compatibility for older files.
        Float64(node["nu_theta_station1"])
    else
        0.0
    end
    m_tip_bounds = haskey(node, "m_tip_bounds") ? (Float64(node["m_tip_bounds"][1]), Float64(node["m_tip_bounds"][2])) : (0.5, 1.2)
    phi_in_bounds = haskey(node, "phi_in_bounds") ? (Float64(node["phi_in_bounds"][1]), Float64(node["phi_in_bounds"][2])) : (0.2, 0.9)
    rows = AxialRow[_row_from_toml_dict(row_data) for row_data in node["rows"]]
    return AxialMachineModel(
        Float64(node["gamma"]),
        Float64(node["gas_constant"]),
        Float64(node["A_ref"]),
        Float64.(node["A_station"]),
        rows,
        nu_theta_first_rotor,
        m_tip_bounds,
        phi_in_bounds,
    )
end

"""
Demo axial-machine model for development/testing.
"""
function demo_axial_machine_model()
    rows = AxialRow[
        AxialRow(
            :rotor,
            RotorAeroModel{Float64}(
                beta_ref=-0.55,
                beta_incidence_sensitivity=0.62,
                loss_base=0.0025,
                loss_incidence=0.045,
                stall_incidence_limit=0.36,
                k_theta_min=-2.0,
                k_theta_max=1.1,
            ),
            0.180,
            0.220,
            +1,
        ),
        AxialRow(
            :stator,
            StatorAeroModel{Float64}(
                alpha_ref=0.45,
                alpha_incidence_sensitivity=0.70,
                loss_base=0.0018,
                loss_incidence=0.030,
                stall_incidence_limit=0.34,
                k_theta_min=-1.0,
                k_theta_max=2.0,
            ),
            0.180,
            0.220,
            0,
        ),
        AxialRow(
            :rotor,
            RotorAeroModel{Float64}(
                beta_ref=-0.50,
                beta_incidence_sensitivity=0.60,
                loss_base=0.0030,
                loss_incidence=0.050,
                stall_incidence_limit=0.34,
                k_theta_min=-2.1,
                k_theta_max=1.1,
            ),
            0.180,
            0.220,
            +1,
        ),
        AxialRow(
            :stator,
            StatorAeroModel{Float64}(
                alpha_ref=0.45,
                alpha_incidence_sensitivity=0.70,
                loss_base=0.0020,
                loss_incidence=0.032,
                stall_incidence_limit=0.34,
                k_theta_min=-1.0,
                k_theta_max=2.0,
            ),
            0.180,
            0.220,
            0,
        ),
    ]
    return AxialMachineModel(
        1.4,
        287.05,
        0.060,
        [1.00, 0.98, 0.96, 0.95, 0.95],
        rows,
        0.0,
        (0.01, 1.10),
        (0.01, 0.95),
    )
end
