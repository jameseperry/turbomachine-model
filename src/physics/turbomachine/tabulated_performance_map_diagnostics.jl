using TOML
import ...Utility: write_toml, read_toml

struct TabulatedPerformanceDiagnosticSample
    omega_corr::Float64
    omega::Float64
    Vx_target::Float64
    Vx_used::Float64
    mdot_corr::Float64
    mdot_target::Float64
    mdot_used::Float64
    exact_feasible::Bool
    repaired::Bool
    result::AxialMachine.StreamtubeSolveResult
end

struct TabulatedPerformanceMapDiagnostics
    Tt_ref::Float64
    Pt_ref::Float64
    speed_grid::Vector{Float64}
    flow_grid::Vector{Float64}
    samples::Matrix{TabulatedPerformanceDiagnosticSample}
end

function TabulatedPerformanceMapDiagnostics(
    Tt_ref::Real,
    Pt_ref::Real,
    speed_grid::AbstractVector{<:Real},
    flow_grid::AbstractVector{<:Real},
    samples::AbstractMatrix{TabulatedPerformanceDiagnosticSample},
)
    size(samples) == (length(speed_grid), length(flow_grid)) ||
        error("samples size must match (length(speed_grid), length(flow_grid))")
    return TabulatedPerformanceMapDiagnostics(
        Float64(Tt_ref),
        Float64(Pt_ref),
        Float64.(speed_grid),
        Float64.(flow_grid),
        Matrix{TabulatedPerformanceDiagnosticSample}(samples),
    )
end

diagnostic_sample(
    diagnostics::TabulatedPerformanceMapDiagnostics,
    idx::CartesianIndex{2},
) = diagnostics.samples[idx]

diagnostic_sample(
    diagnostics::TabulatedPerformanceMapDiagnostics,
    i_speed::Integer,
    i_flow::Integer,
) = diagnostics.samples[i_speed, i_flow]

function nearest_diagnostic_samples(
    map::TabulatedPerformanceMap,
    diagnostics::TabulatedPerformanceMapDiagnostics,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real;
    n::Int=5,
)
    diagnostics.speed_grid == tabulated_speed_grid(map) || error("diagnostics speed grid must match map speed grid")
    diagnostics.flow_grid == tabulated_flow_grid(map) || error("diagnostics flow grid must match map flow grid")
    indices = nearest_grid_indices(map, omega, mdot, Tt_in, Pt_in; n=n)
    return [diagnostics.samples[idx] for idx in indices]
end

function _station_state_to_toml_dict(st::AxialMachine.StreamtubeStationState)
    return Dict{String,Any}(
        "radius" => st.radius,
        "area" => st.area,
        "pt_t" => st.pt_t,
        "ht_t" => st.ht_t,
        "s_t" => st.s_t,
        "p" => st.p,
        "h" => st.h,
        "Tt" => st.Tt,
        "T" => st.T,
        "rho" => st.rho,
        "a" => st.a,
        "Vx" => st.Vx,
        "Vtheta" => st.Vtheta,
        "V" => st.V,
        "Mach" => st.Mach,
        "mdot_station" => st.mdot_station,
        "valid" => st.valid,
        "stall" => st.stall,
        "choke" => st.choke,
    )
end

function _station_state_from_toml_dict(data::Dict{String,Any})
    return AxialMachine.StreamtubeStationState(
        Float64(data["radius"]),
        Float64(data["area"]),
        Float64(data["pt_t"]),
        Float64(data["ht_t"]),
        Float64(data["s_t"]),
        Float64(data["p"]),
        Float64(data["h"]),
        Float64(data["Tt"]),
        Float64(data["T"]),
        Float64(data["rho"]),
        Float64(data["a"]),
        Float64(data["Vx"]),
        Float64(data["Vtheta"]),
        Float64(data["V"]),
        Float64(data["Mach"]),
        Float64(data["mdot_station"]),
        Bool(data["valid"]),
        Bool(data["stall"]),
        Bool(data["choke"]),
    )
end

function _aero_result_to_toml_dict(aero::AxialMachine.BladeAeroResult)
    return Dict{String,Any}(
        "theta_out" => aero.theta_out,
        "theta_out_unclamped" => aero.theta_out_unclamped,
        "delta_s_t" => aero.delta_s_t,
        "stall_margin" => aero.stall_margin,
        "regime" => String(aero.regime),
        "incidence" => aero.incidence,
        "deviation" => aero.deviation,
        "theta_in" => aero.theta_in,
        "theta_metal_in" => aero.theta_metal_in,
        "theta_metal_out" => aero.theta_metal_out,
    )
end

function _aero_result_from_toml_dict(data::Dict{String,Any})
    return AxialMachine.BladeAeroResult(
        Float64(data["theta_out"]),
        Float64(data["theta_out_unclamped"]),
        Float64(data["delta_s_t"]),
        Float64(data["stall_margin"]),
        Symbol(String(data["regime"])),
        Float64(data["incidence"]),
        Float64(data["deviation"]),
        Float64(data["theta_in"]),
        Float64(data["theta_metal_in"]),
        Float64(data["theta_metal_out"]),
    )
end

function _row_state_to_toml_dict(row::AxialMachine.StreamtubeRowState)
    return Dict{String,Any}(
        "row_index" => row.row_index,
        "aero" => _aero_result_to_toml_dict(row.aero),
        "selected_candidate_index" => row.selected_candidate_index,
        "choke" => row.choke,
        "outlet_candidates" => [_station_state_to_toml_dict(st) for st in row.outlet_candidates],
    )
end

function _row_state_from_toml_dict(data::Dict{String,Any})
    return AxialMachine.StreamtubeRowState(
        row_index=Int(data["row_index"]),
        aero=_aero_result_from_toml_dict(data["aero"]),
        outlet_candidates=[_station_state_from_toml_dict(st) for st in data["outlet_candidates"]],
        selected_candidate_index=Int(data["selected_candidate_index"]),
        choke=Bool(data["choke"]),
    )
end

function _streamtube_result_to_toml_dict(result::AxialMachine.StreamtubeSolveResult)
    return Dict{String,Any}(
        "pressure_ratio" => result.pressure_ratio,
        "efficiency" => result.efficiency,
        "thermo_efficiency" => result.thermo_efficiency,
        "valid" => result.valid,
        "stall" => result.stall,
        "choke" => result.choke,
        "status" => String(result.status),
        "stations" => [_station_state_to_toml_dict(st) for st in result.stations],
        "row_data" => [_row_state_to_toml_dict(row) for row in result.row_data],
    )
end

function _streamtube_result_from_toml_dict(data::Dict{String,Any})
    return AxialMachine.StreamtubeSolveResult(
        Float64(data["pressure_ratio"]),
        Float64(data["efficiency"]),
        Float64(data["thermo_efficiency"]),
        Bool(data["valid"]),
        Bool(data["stall"]),
        Bool(data["choke"]),
        Symbol(String(data["status"])),
        [_station_state_from_toml_dict(st) for st in data["stations"]],
        [_row_state_from_toml_dict(row) for row in data["row_data"]],
    )
end

function _diagnostic_sample_to_toml_dict(sample::TabulatedPerformanceDiagnosticSample, i_speed::Int, i_flow::Int)
    return Dict{String,Any}(
        "i_speed" => i_speed,
        "i_flow" => i_flow,
        "omega_corr" => sample.omega_corr,
        "omega" => sample.omega,
        "Vx_target" => sample.Vx_target,
        "Vx_used" => sample.Vx_used,
        "mdot_corr" => sample.mdot_corr,
        "mdot_target" => sample.mdot_target,
        "mdot_used" => sample.mdot_used,
        "exact_feasible" => sample.exact_feasible,
        "repaired" => sample.repaired,
        "result" => _streamtube_result_to_toml_dict(sample.result),
    )
end

function _diagnostic_sample_from_toml_dict(data::Dict{String,Any})
    return (
        i_speed=Int(data["i_speed"]),
        i_flow=Int(data["i_flow"]),
        sample=TabulatedPerformanceDiagnosticSample(
            Float64(data["omega_corr"]),
            Float64(data["omega"]),
            Float64(data["Vx_target"]),
            Float64(data["Vx_used"]),
            Float64(data["mdot_corr"]),
            Float64(data["mdot_target"]),
            Float64(data["mdot_used"]),
            Bool(data["exact_feasible"]),
            Bool(data["repaired"]),
            _streamtube_result_from_toml_dict(data["result"]),
        ),
    )
end

function write_toml(
    diagnostics::TabulatedPerformanceMapDiagnostics,
    path::AbstractString;
    group::AbstractString="performance_map_diagnostics",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "tabulated_performance_map_diagnostics"
    node["format_version"] = 1
    node["Tt_ref"] = diagnostics.Tt_ref
    node["Pt_ref"] = diagnostics.Pt_ref
    node["speed_grid"] = diagnostics.speed_grid
    node["flow_grid"] = diagnostics.flow_grid
    node["samples"] = vec([
        _diagnostic_sample_to_toml_dict(diagnostics.samples[i, j], i, j)
        for i in axes(diagnostics.samples, 1), j in axes(diagnostics.samples, 2)
    ])
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{TabulatedPerformanceMapDiagnostics},
    path::AbstractString;
    group::AbstractString="performance_map_diagnostics",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    speed_grid = Float64.(node["speed_grid"])
    flow_grid = Float64.(node["flow_grid"])
    samples = Matrix{TabulatedPerformanceDiagnosticSample}(undef, length(speed_grid), length(flow_grid))
    for sample_data in node["samples"]
        parsed = _diagnostic_sample_from_toml_dict(sample_data)
        samples[parsed.i_speed, parsed.i_flow] = parsed.sample
    end
    return TabulatedPerformanceMapDiagnostics(
        Float64(node["Tt_ref"]),
        Float64(node["Pt_ref"]),
        speed_grid,
        flow_grid,
        samples,
    )
end
