#!/usr/bin/env julia

using ArgParse
using Statistics
using TurboMachineModel

const TM = TurboMachineModel.Physics.Turbomachine
const Axial = TM.AxialMachine

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _pctl(values::Vector{Float64}, p::Float64)
    isempty(values) && return NaN
    s = sort(values)
    idx = clamp(Int(ceil(p * length(s))), 1, length(s))
    return s[idx]
end

function _finite_values(values::Vector{Float64})
    return [v for v in values if isfinite(v)]
end

function _roughness_metric(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || error("x/y length mismatch")
    n >= 3 || return NaN
    all(isfinite, y) || return NaN

    slopes = Float64[]
    for i in 1:(n - 1)
        dx = x[i + 1] - x[i]
        dx > 0 || continue
        push!(slopes, (y[i + 1] - y[i]) / dx)
    end
    length(slopes) >= 2 || return NaN

    curv = Float64[]
    for i in 1:(length(slopes) - 1)
        push!(curv, abs(slopes[i + 1] - slopes[i]))
    end
    isempty(curv) && return NaN

    y_rng = max(maximum(y) - minimum(y), 1e-12)
    return median(curv) / y_rng
end

function _stats_line(label::String, vals::Vector{Float64})
    f = _finite_values(vals)
    if isempty(f)
        println("  $(label): no finite samples")
        return
    end
    println("  $(label): n=$(length(f)), median=$(median(f)), p95=$(_pctl(f, 0.95)), max=$(maximum(f))")
end

function _eta_turb_from_comp(eta_comp::Float64)
    if eta_comp < -1e-8
        return clamp(-1.0 / eta_comp, 0.0, 1.0)
    end
    return NaN
end

function _direct_turbine_state_from_mdot(
    model::Axial.AxialMachineModel,
    m_tip::Float64,
    mdot::Float64,
    Tt_in_ref::Float64,
    Pt_in_ref::Float64,
)
    eos = TurboMachineModel.Physics.Fluids.IdealGasEOS(:axial; gas_constant=model.gas_constant, gamma=model.gamma)
    ht_in_ref = TurboMachineModel.Physics.Fluids.enthalpy_from_temperature(eos, Tt_in_ref)
    a0_in_ref = sqrt(model.gamma * model.gas_constant * Tt_in_ref)
    omega = m_tip * a0_in_ref / model.r_tip_ref
    vals = Axial.streamtube_solve_from_mdot(
        model,
        eos,
        Axial.meanline_radii(model),
        Pt_in_ref,
        ht_in_ref,
        omega,
        mdot,
        0.0,
    )
    vals.valid || return (valid=false,)
    isfinite(vals.PR) || return (valid=false,)
    vals.PR > 0 || return (valid=false,)
    vals.PR < 1 || return (valid=false,)

    eta_t = Float64(vals.eta)
    isfinite(eta_t) || return (valid=false,)

    return (
        valid=true,
        PR=vals.PR,
        eta=eta_t,
        omega=omega,
        mdot=mdot,
    )
end

function _build_parser()
    settings = ArgParseSettings(
        prog="diagnose_turbine_map_pipeline.jl",
        description="Diagnose whether turbine map roughness originates in direct streamtube samples or common-map tabulation.",
    )
    @add_arg_table! settings begin
        "meanline_path"
            help = "input axial meanline model TOML"
            required = true
        "map_path"
            help = "input common performance map TOML"
            required = true
        "--meanline-group"
            help = "input meanline TOML group"
            arg_type = String
            default = "axial_model"
        "--map-group"
            help = "input map TOML group"
            arg_type = String
            default = "performance_map"
        "--tt-in-ref"
            help = "reference inlet total temperature (K)"
            arg_type = Float64
            default = 288.15
        "--pt-in-ref"
            help = "reference inlet total pressure (Pa)"
            arg_type = Float64
            default = 101_325.0
        "--n-speed"
            help = "diagnostic speed samples in m_tip space"
            arg_type = Int
            default = 31
        "--n-phi"
            help = "diagnostic flow samples in mdot space"
            arg_type = Int
            default = 121
        "--csv"
            help = "optional per-sample comparison CSV path"
            arg_type = String
    end
    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())

    meanline_group = String(something(_parsed_opt(parsed, "meanline_group", "meanline-group"), "axial_model"))
    map_group = String(something(_parsed_opt(parsed, "map_group", "map-group"), "performance_map"))
    Tt_in_ref = Float64(something(_parsed_opt(parsed, "tt_in_ref", "tt-in-ref"), 288.15))
    Pt_in_ref = Float64(something(_parsed_opt(parsed, "pt_in_ref", "pt-in-ref"), 101_325.0))
    n_speed = Int(something(_parsed_opt(parsed, "n_speed", "n-speed"), 31))
    n_phi = Int(something(_parsed_opt(parsed, "n_phi", "n-phi"), 121))
    csv_path = _parsed_opt(parsed, "csv", "csv")

    model = Axial.read_toml(Axial.AxialMachineModel, parsed["meanline_path"]; group=meanline_group)
    map = TM.read_toml(TM.TabulatedPerformanceMap, parsed["map_path"]; group=map_group)

    m_lo, m_hi = model.m_tip_bounds
    m_grid = collect(range(m_lo, m_hi, length=n_speed))
    eos = TurboMachineModel.Physics.Fluids.IdealGasEOS(:axial; gas_constant=model.gas_constant, gamma=model.gamma)
    ht_in_ref = TurboMachineModel.Physics.Fluids.enthalpy_from_temperature(eos, Tt_in_ref)
    inlet_area = Axial.station_area(model, 1)
    mdot_grid = let
        Vx_grid = collect(range(model.Vx_bounds[1], model.Vx_bounds[2], length=n_phi))
        [
            Axial._build_station_state(
                eos,
                Axial._station_radius(model, Axial.meanline_radii(model), 1),
                inlet_area,
                Pt_in_ref,
                ht_in_ref,
                Vx,
                0.0,
                NaN,
                true,
                false,
                false,
            ).mdot_station for Vx in Vx_grid
        ]
    end

    pr_err = Float64[]
    eta_err = Float64[]
    direct_pr_rough = Float64[]
    direct_eta_rough = Float64[]
    map_pr_rough = Float64[]
    map_eta_rough = Float64[]
    sample_rows = NamedTuple[]

    for (i, m_tip) in pairs(m_grid)
        direct_line = NamedTuple[]
        map_line = NamedTuple[]
        for mdot in mdot_grid
            st = _direct_turbine_state_from_mdot(model, Float64(m_tip), Float64(mdot), Tt_in_ref, Pt_in_ref)
            st.valid || continue
            map_vals = TM.performance_from_stagnation(map, st.omega, st.mdot, Tt_in_ref, Pt_in_ref)
            map_vals.valid || continue

            row = (
                speed_idx=i,
                m_tip=Float64(m_tip),
                mdot_input=Float64(mdot),
                omega=st.omega,
                mdot=st.mdot,
                pr_direct=st.PR,
                eta_direct=st.eta,
                pr_map=Float64(map_vals.pressure_ratio),
                eta_map=Float64(map_vals.eta),
                pr_err=abs(Float64(map_vals.pressure_ratio) - st.PR),
                eta_err=abs(Float64(map_vals.eta) - st.eta),
            )
            push!(sample_rows, row)
            push!(direct_line, (x=st.mdot, pr=st.PR, eta=st.eta))
            push!(map_line, (x=st.mdot, pr=Float64(map_vals.pressure_ratio), eta=Float64(map_vals.eta)))
            push!(pr_err, row.pr_err)
            push!(eta_err, row.eta_err)
        end

        if length(direct_line) >= 3
            direct_sorted = sort(direct_line; by=r -> r.x)
            map_sorted = sort(map_line; by=r -> r.x)
            x_d = [r.x for r in direct_sorted]
            pr_d = [r.pr for r in direct_sorted]
            eta_d = [r.eta for r in direct_sorted]
            x_m = [r.x for r in map_sorted]
            pr_m = [r.pr for r in map_sorted]
            eta_m = [r.eta for r in map_sorted]
            push!(direct_pr_rough, _roughness_metric(x_d, pr_d))
            push!(direct_eta_rough, _roughness_metric(x_d, eta_d))
            push!(map_pr_rough, _roughness_metric(x_m, pr_m))
            push!(map_eta_rough, _roughness_metric(x_m, eta_m))
        end
    end

    println("Turbine Map Pipeline Diagnostics")
    println("- meanline_path: $(parsed["meanline_path"])")
    println("- map_path: $(parsed["map_path"])")
    _stats_line("PR abs error", pr_err)
    _stats_line("eta abs error", eta_err)
    _stats_line("direct PR roughness", direct_pr_rough)
    _stats_line("map PR roughness", map_pr_rough)
    _stats_line("direct eta roughness", direct_eta_rough)
    _stats_line("map eta roughness", map_eta_rough)

    if !isnothing(csv_path)
        open(csv_path, "w") do io
            println(io, "speed_idx,m_tip,mdot_input,omega,mdot,pr_direct,eta_direct,pr_map,eta_map,pr_err,eta_err")
            for row in sample_rows
                println(io, join((
                    row.speed_idx,
                    row.m_tip,
                    row.mdot_input,
                    row.omega,
                    row.mdot,
                    row.pr_direct,
                    row.eta_direct,
                    row.pr_map,
                    row.eta_map,
                    row.pr_err,
                    row.eta_err,
                ), ','))
            end
        end
        println("Wrote diagnostic CSV to: $(csv_path)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
