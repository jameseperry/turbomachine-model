using TOML
import ....Utility: write_toml, read_toml

function _aero_to_toml_dict(model::BladeAeroModel)
    return Dict{String,Any}(
        "format" => "blade_aero_model",
        "theta_ref" => Float64(model.theta_ref),
        "theta_incidence_sensitivity" => Float64(model.theta_incidence_sensitivity),
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
    if fmt == "blade_aero_model"
        return BladeAeroModel{Float64}(
            Float64(data["theta_ref"]),
            Float64(data["theta_incidence_sensitivity"]),
            Float64(data["loss_base"]),
            Float64(data["loss_incidence"]),
            Float64(data["stall_incidence_limit"]),
            Float64(data["k_theta_min"]),
            Float64(data["k_theta_max"]),
        )
    end
    error("unsupported aero model format=$(fmt)")
end

function _row_to_toml_dict(row::AxialRow)
    return Dict{String,Any}(
        "r_hub" => row.r_hub,
        "r_tip" => row.r_tip,
        "speed_ratio_to_ref" => row.speed_ratio_to_ref,
        "aero" => _aero_to_toml_dict(row.aero),
    )
end

function _row_from_toml_dict(data::Dict{String,Any})
    haskey(data, "r_hub") || error("row missing r_hub")
    haskey(data, "r_tip") || error("row missing r_tip")
    haskey(data, "speed_ratio_to_ref") || error("row missing speed_ratio_to_ref")
    haskey(data, "aero") || error("row missing aero")
    aero = _aero_from_toml_dict(data["aero"])
    return AxialRow(
        aero,
        Float64(data["r_hub"]),
        Float64(data["r_tip"]),
        Float64(data["speed_ratio_to_ref"]),
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
    node["format"] = "compressor_meanline_model"
    node["format_version"] = 5
    node["gamma"] = model.gamma
    node["gas_constant"] = model.gas_constant
    node["r_tip_ref"] = model.r_tip_ref
    node["r_flow_ref"] = model.r_flow_ref
    node["speed_ratio_ref"] = model.speed_ratio_ref
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
    for key in ("gamma", "gas_constant", "r_tip_ref", "r_flow_ref", "speed_ratio_ref", "m_tip_bounds", "phi_in_bounds", "rows")
        haskey(node, key) || error("missing TOML key $(key)")
    end
    rows = AxialRow[_row_from_toml_dict(row_data) for row_data in node["rows"]]
    return AxialMachineModel(
        Float64(node["gamma"]),
        Float64(node["gas_constant"]),
        Float64(node["r_tip_ref"]),
        rows,
        (Float64(node["m_tip_bounds"][1]), Float64(node["m_tip_bounds"][2])),
        (Float64(node["phi_in_bounds"][1]), Float64(node["phi_in_bounds"][2])),
        r_flow_ref=Float64(node["r_flow_ref"]),
        speed_ratio_ref=Float64(node["speed_ratio_ref"]),
    )
end

"""
Demo axial-machine model for development/testing.
"""
function demo_axial_compressor_model()
    rows = AxialRow[
        AxialRow(
            rotor_aero_model(
                theta_ref=-0.55,
                theta_incidence_sensitivity=0.62,
                loss_base=0.0025,
                loss_incidence=0.045,
                stall_incidence_limit=0.36,
                k_theta_min=-2.0,
                k_theta_max=1.1,
            ),
            0.140,
            0.220,
            1.0,
        ),
        AxialRow(
            stator_aero_model(
                theta_ref=0.45,
                theta_incidence_sensitivity=0.70,
                loss_base=0.0018,
                loss_incidence=0.030,
                stall_incidence_limit=0.34,
                k_theta_min=-1.0,
                k_theta_max=2.0,
            ),
            0.140,
            0.220,
            0.0,
        ),
        AxialRow(
            rotor_aero_model(
                theta_ref=-0.50,
                theta_incidence_sensitivity=0.60,
                loss_base=0.0030,
                loss_incidence=0.050,
                stall_incidence_limit=0.34,
                k_theta_min=-2.1,
                k_theta_max=1.1,
            ),
            0.140,
            0.220,
            1.0,
        ),
        AxialRow(
            stator_aero_model(
                theta_ref=0.45,
                theta_incidence_sensitivity=0.70,
                loss_base=0.0020,
                loss_incidence=0.032,
                stall_incidence_limit=0.34,
                k_theta_min=-1.0,
                k_theta_max=2.0,
            ),
            0.140,
            0.220,
            0.0,
        ),
    ]
    return AxialMachineModel(
        1.4,
        287.05,
        0.220,
        rows,
        (0.01, 1.10),
        (0.01, 0.95),
    )
end

"""
Demo axial-machine model configured to behave turbine-like for development/testing.

This model is still serialized using the same `compressor_meanline_model` schema as
the compressor demo; only the aero row parameters differ.
"""
function demo_axial_turbine_model()
    rows = AxialRow[
        # Simplified single-stage turbine-like setup: one guide vane and one rotor.
        AxialRow(
            stator_aero_model(
                theta_ref=-0.55,
                theta_incidence_sensitivity=0.35,
                loss_base=0.0010,
                loss_incidence=0.010,
                stall_incidence_limit=0.55,
                k_theta_min=-1.8,
                k_theta_max=0.5,
            ),
            0.140,
            0.220,
            0.0,
        ),
        # Rotor row tuned for smoother work extraction over a narrower domain.
        AxialRow(
            rotor_aero_model(
                theta_ref=0.60,
                theta_incidence_sensitivity=0.35,
                loss_base=0.0012,
                loss_incidence=0.012,
                stall_incidence_limit=0.55,
                k_theta_min=-0.5,
                k_theta_max=1.8,
            ),
            0.140,
            0.220,
            1.0,
        ),
    ]
    return AxialMachineModel(
        1.4,
        287.05,
        0.220,
        rows,
        (0.10, 0.90),
        (0.10, 0.70),
    )
end
