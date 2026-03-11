using TOML
import ....Utility: write_toml, read_toml

function _aero_to_toml_dict(model::BladeAeroModel)
    return Dict{String,Any}(
        "format" => "blade_aero_model",
        "deviation_ref" => Float64(model.deviation_ref),
        "deviation_incidence_sensitivity" => Float64(model.deviation_incidence_sensitivity),
        "loss_entropy_base" => Float64(model.loss_entropy_base),
        "loss_entropy_incidence" => Float64(model.loss_entropy_incidence),
        "stall_incidence_limit" => Float64(model.stall_incidence_limit),
        "theta_min" => Float64(model.theta_min),
        "theta_max" => Float64(model.theta_max),
    )
end

function _aero_from_toml_dict(data::Dict{String,Any})
    haskey(data, "format") || error("aero model missing format")
    fmt = String(data["format"])
    if fmt == "blade_aero_model"
        return BladeAeroModel{Float64}(
            Float64(data["deviation_ref"]),
            Float64(data["deviation_incidence_sensitivity"]),
            Float64(data["loss_entropy_base"]),
            Float64(data["loss_entropy_incidence"]),
            Float64(data["stall_incidence_limit"]),
            Float64(data["theta_min"]),
            Float64(data["theta_max"]),
        )
    end
    error("unsupported aero model format=$(fmt)")
end

function _row_to_toml_dict(row::AxialRow)
    return Dict{String,Any}(
        "r_hub" => row.r_hub,
        "r_tip" => row.r_tip,
        "theta_metal_in" => row.theta_metal_in,
        "theta_metal_out" => row.theta_metal_out,
        "speed_ratio_to_ref" => row.speed_ratio_to_ref,
        "aero" => _aero_to_toml_dict(row.aero),
    )
end

function _row_from_toml_dict(data::Dict{String,Any})
    haskey(data, "r_hub") || error("row missing r_hub")
    haskey(data, "r_tip") || error("row missing r_tip")
    haskey(data, "theta_metal_in") || error("row missing theta_metal_in")
    haskey(data, "theta_metal_out") || error("row missing theta_metal_out")
    haskey(data, "speed_ratio_to_ref") || error("row missing speed_ratio_to_ref")
    haskey(data, "aero") || error("row missing aero")
    aero = _aero_from_toml_dict(data["aero"])
    return AxialRow(
        aero,
        Float64(data["r_hub"]),
        Float64(data["r_tip"]),
        Float64(data["theta_metal_in"]),
        Float64(data["theta_metal_out"]),
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
    group::AbstractString="axial_model",
)
    data = Dict{String,Any}()
    node = _find_or_create_axial_group!(data, group)
    node["format"] = "axial_model"
    node["format_version"] = 6
    node["gamma"] = model.gamma
    node["gas_constant"] = model.gas_constant
    node["r_tip_ref"] = model.r_tip_ref
    node["r_flow_ref"] = model.r_flow_ref
    node["speed_ratio_ref"] = model.speed_ratio_ref
    node["m_tip_bounds"] = [model.m_tip_bounds[1], model.m_tip_bounds[2]]
    node["Vx_bounds"] = [model.Vx_bounds[1], model.Vx_bounds[2]]
    node["rows"] = [_row_to_toml_dict(row) for row in model.rows]
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{AxialMachineModel},
    path::AbstractString;
    group::AbstractString="axial_model",
)
    data = TOML.parsefile(path)
    node = _find_axial_group(data, group)
    for key in ("gamma", "gas_constant", "r_tip_ref", "r_flow_ref", "speed_ratio_ref", "m_tip_bounds", "Vx_bounds", "rows")
        haskey(node, key) || error("missing TOML key $(key)")
    end
    rows = AxialRow[_row_from_toml_dict(row_data) for row_data in node["rows"]]
    return AxialMachineModel(
        Float64(node["gamma"]),
        Float64(node["gas_constant"]),
        Float64(node["r_tip_ref"]),
        rows,
        (Float64(node["m_tip_bounds"][1]), Float64(node["m_tip_bounds"][2])),
        (Float64(node["Vx_bounds"][1]), Float64(node["Vx_bounds"][2])),
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
                deviation_ref=0.0,
                deviation_incidence_sensitivity=0.62,
                loss_entropy_base=2.5,
                loss_entropy_incidence=45.0,
                stall_incidence_limit=0.36,
                theta_min=atan(-2.0),
                theta_max=atan(1.1),
            ),
            0.140,
            0.220,
            -0.55,
            -0.55,
            1.0,
        ),
        AxialRow(
            stator_aero_model(
                deviation_ref=0.0,
                deviation_incidence_sensitivity=0.70,
                loss_entropy_base=1.8,
                loss_entropy_incidence=30.0,
                stall_incidence_limit=0.34,
                theta_min=atan(-1.0),
                theta_max=atan(2.0),
            ),
            0.140,
            0.220,
            0.45,
            0.45,
            0.0,
        ),
        AxialRow(
            rotor_aero_model(
                deviation_ref=0.0,
                deviation_incidence_sensitivity=0.45,
                loss_entropy_base=2.4,
                loss_entropy_incidence=36.0,
                stall_incidence_limit=0.34,
                theta_min=atan(-1.5),
                theta_max=atan(1.1),
            ),
            0.140,
            0.220,
            -0.40,
            -0.30,
            1.0,
        ),
        AxialRow(
            stator_aero_model(
                deviation_ref=0.0,
                deviation_incidence_sensitivity=0.70,
                loss_entropy_base=2.0,
                loss_entropy_incidence=32.0,
                stall_incidence_limit=0.34,
                theta_min=atan(-1.0),
                theta_max=atan(2.0),
            ),
            0.140,
            0.220,
            0.45,
            0.45,
            0.0,
        ),
    ]
    return AxialMachineModel(
        1.4,
        287.05,
        0.220,
        rows,
        (0.01, 1.10),
        (4.0, 220.0),
    )
end

"""
Demo axial-machine model configured to behave turbine-like for development/testing.

This model is serialized using the common `axial_model` schema; only the aero
row parameters differ from the compressor demo.
"""
function demo_axial_turbine_model()
    rows = AxialRow[
        # Simplified single-stage turbine-like setup: one guide vane and one rotor.
        AxialRow(
            stator_aero_model(
                deviation_ref=0.0,
                deviation_incidence_sensitivity=0.45,
                loss_entropy_base=1.0,
                loss_entropy_incidence=10.0,
                stall_incidence_limit=0.55,
                theta_min=atan(0.1),
                theta_max=atan(2.2),
            ),
            0.140,
            0.220,
            0.0,
            deg2rad(70.0),
            0.0,
        ),
        # Rotor row tuned for smoother work extraction over a narrower domain.
        AxialRow(
            rotor_aero_model(
                deviation_ref=0.0,
                deviation_incidence_sensitivity=0.55,
                loss_entropy_base=1.2,
                loss_entropy_incidence=12.0,
                stall_incidence_limit=0.55,
                theta_min=atan(-3.0),
                theta_max=atan(-0.6),
            ),
            0.140,
            0.220,
            deg2rad(40.0),
            deg2rad(-50.0),
            1.0,
        ),
    ]
    return AxialMachineModel(
        1.4,
        287.05,
        0.220,
        rows,
        (0.10, 0.90),
        (6.0, 120.0),
    )
end
