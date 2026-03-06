#!/usr/bin/env julia

using ArgParse
using Statistics
using TurboMachineModel

const Axial = TurboMachineModel.Physics.Turbomachine.AxialMachine
const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine

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
    y_finite = [isfinite(v) for v in y]
    all(y_finite) || return NaN

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
    println(
        "  $(label): n=$(length(f)), median=$(median(f)), p95=$(_pctl(f, 0.95)), max=$(maximum(f))",
    )
end

function _eta_turb_from_comp(eta_comp::Float64)
    if eta_comp < -1e-8
        return clamp(-1.0 / eta_comp, 0.0, 1.0)
    end
    return NaN
end

function _direct_turbine_state_from_phi(
    model::Axial.AxialMachineModel,
    m_tip::Float64,
    phi_in::Float64,
    Tt_in_ref::Float64,
    Pt_in_ref::Float64,
    Tt_ref::Float64,
    Pt_ref::Float64,
)
    vals = Axial.streamtube_solve_with_phi(model, m_tip, phi_in)
    vals.valid || return (valid=false,)
    isfinite(vals.PR) || return (valid=false,)
    vals.PR > 0 || return (valid=false,)
    vals.PR < 1 || return (valid=false,)

    eta_t = _eta_turb_from_comp(Float64(vals.eta))
    isfinite(eta_t) || return (valid=false,)

    a0_in_ref = sqrt(model.gamma * model.gas_constant * Tt_in_ref)
    omega = m_tip * a0_in_ref / model.r_tip_ref
    omega_corr = Turbine.corrected_speed(omega, Tt_in_ref, Tt_ref)
    rho0_in_ref = Pt_in_ref / (model.gas_constant * Tt_in_ref)
    inlet_area = Axial.station_area(model, 1)
    mean_radius_inlet = abs(model.speed_ratio_ref) * model.r_flow_ref
    mdot = phi_in * rho0_in_ref * inlet_area * abs(omega * mean_radius_inlet)
    mdot_corr = Turbine.corrected_flow(mdot, Tt_in_ref, Pt_in_ref, Tt_ref, Pt_ref)

    return (
        valid=true,
        PR_turb=1.0 / vals.PR,
        eta=eta_t,
        mdot_corr=mdot_corr,
        omega_corr=omega_corr,
    )
end

function _build_parser()
    settings = ArgParseSettings(
        prog="diagnose_turbine_map_pipeline.jl",
        description="Diagnose whether turbine map roughness originates in streamtube samples, tabulation, or interpolation.",
    )
    @add_arg_table! settings begin
        "meanline_path"
            help = "input axial meanline model TOML"
            required = true
        "map_path"
            help = "input tabulated turbine map TOML"
            required = true
        "--meanline-group"
            help = "input meanline TOML group"
            arg_type = String
            default = "compressor_meanline_model"
        "--map-group"
            help = "input map TOML group"
            arg_type = String
            default = "turbine_map"
        "--tt-in-ref"
            help = "reference inlet total temperature used for meanline-to-map conversion (K)"
            arg_type = Float64
            default = 288.15
        "--pt-in-ref"
            help = "reference inlet total pressure used for meanline-to-map conversion (Pa)"
            arg_type = Float64
            default = 101_325.0
        "--n-speed"
            help = "diagnostic speed samples in m_tip space"
            arg_type = Int
            default = 31
        "--n-phi"
            help = "diagnostic flow samples in phi space"
            arg_type = Int
            default = 121
        "--n-pr-interp"
            help = "interpolated PR samples per speed for interpolation roughness checks"
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

    meanline_group = something(_parsed_opt(parsed, "meanline_group", "meanline-group"), "compressor_meanline_model")
    map_group = something(_parsed_opt(parsed, "map_group", "map-group"), "turbine_map")
    Tt_in_ref = Float64(something(_parsed_opt(parsed, "tt_in_ref", "tt-in-ref"), 288.15))
    Pt_in_ref = Float64(something(_parsed_opt(parsed, "pt_in_ref", "pt-in-ref"), 101_325.0))
    n_speed = Int(something(_parsed_opt(parsed, "n_speed", "n-speed"), 31))
    n_phi = Int(something(_parsed_opt(parsed, "n_phi", "n-phi"), 121))
    n_pr_interp = Int(something(_parsed_opt(parsed, "n_pr_interp", "n-pr-interp"), 121))
    csv_path = _parsed_opt(parsed, "csv", "csv")

    model = Axial.read_toml(Axial.AxialMachineModel, parsed["meanline_path"]; group=meanline_group)
    map = Turbine.read_toml(Turbine.TabulatedTurbinePerformanceMap, parsed["map_path"]; group=map_group)

    m_lo, m_hi = model.m_tip_bounds
    phi_lo, phi_hi = model.phi_in_bounds
    m_grid = collect(range(m_lo, m_hi, length=n_speed))
    phi_grid = collect(range(phi_lo, phi_hi, length=n_phi))

    Tt_ref = map.Tt_ref
    Pt_ref = map.Pt_ref

    direct_mdot_rough = Float64[]
    direct_eta_rough = Float64[]
    table_mdot_rough = Float64[]
    table_eta_rough = Float64[]
    interp_mdot_rough = Float64[]
    interp_eta_rough = Float64[]
    direct_eta_rough_by_speed = NamedTuple[]
    table_eta_rough_by_speed = NamedTuple[]
    mdot_err = Float64[]
    eta_err = Float64[]

    sample_rows = NamedTuple[]

    omega_grid = Turbine.omega_corr_grid(map)
    pr_grid = Turbine.pr_turb_grid(map)
    mdot_table = Turbine.mdot_corr_table(map)
    eta_table = Turbine.eta_table(map)

    for (i, m_tip) in pairs(m_grid)
        direct_line = NamedTuple[]
        for phi in phi_grid
            st = _direct_turbine_state_from_phi(
                model,
                Float64(m_tip),
                Float64(phi),
                Tt_in_ref,
                Pt_in_ref,
                Tt_ref,
                Pt_ref,
            )
            st.valid || continue
            map_vals = Turbine.turbine_performance_map(map, st.omega_corr, st.PR_turb)
            map_valid = hasproperty(map_vals, :valid) ? map_vals.valid : true
            mdot_e = map_valid ? abs(Float64(map_vals.mdot_corr) - st.mdot_corr) : NaN
            eta_e = map_valid ? abs(Float64(map_vals.eta) - st.eta) : NaN

            if isfinite(mdot_e)
                push!(mdot_err, mdot_e)
            end
            if isfinite(eta_e)
                push!(eta_err, eta_e)
            end

            row = (
                speed_idx=i,
                m_tip=Float64(m_tip),
                phi=Float64(phi),
                omega_corr=st.omega_corr,
                PR_turb=st.PR_turb,
                mdot_corr_direct=st.mdot_corr,
                eta_direct=st.eta,
                map_valid=map_valid,
                mdot_corr_map=Float64(map_vals.mdot_corr),
                eta_map=Float64(map_vals.eta),
                mdot_err=mdot_e,
                eta_err=eta_e,
            )
            push!(sample_rows, row)
            push!(direct_line, row)
        end

        if length(direct_line) >= 3
            sorted = sort(direct_line; by=r -> r.PR_turb)
            x = [r.PR_turb for r in sorted]
            md = [r.mdot_corr_direct for r in sorted]
            et = [r.eta_direct for r in sorted]
            r_md = _roughness_metric(x, md)
            r_et = _roughness_metric(x, et)
            push!(direct_mdot_rough, r_md)
            push!(direct_eta_rough, r_et)
            push!(direct_eta_rough_by_speed, (speed_idx=i, m_tip=Float64(m_tip), rough_eta=r_et))
        end
    end

    for i in eachindex(omega_grid)
        lo = map.pr_turb_min[i]
        hi = map.pr_turb_max[i]
        lo > hi && ((lo, hi) = (hi, lo))

        mask = [pr >= lo && pr <= hi for pr in pr_grid]
        idx = findall(identity, mask)
        if length(idx) >= 3
            x = Float64.(pr_grid[idx])
            md = Float64.(mdot_table[i, idx])
            et = Float64.(eta_table[i, idx])
            r_md = _roughness_metric(x, md)
            r_et = _roughness_metric(x, et)
            push!(table_mdot_rough, r_md)
            push!(table_eta_rough, r_et)
            push!(table_eta_rough_by_speed, (speed_idx=i, omega_corr=Float64(omega_grid[i]), rough_eta=r_et))
        end

        pr_dense = collect(range(lo, hi, length=n_pr_interp))
        md_dense = Float64[]
        et_dense = Float64[]
        for pr in pr_dense
            vals = Turbine.turbine_performance_map(map, omega_grid[i], pr)
            if hasproperty(vals, :valid) && !vals.valid
                continue
            end
            push!(md_dense, Float64(vals.mdot_corr))
            eta_val = Float64(vals.eta)
            push!(et_dense, eta_val > 0 && isfinite(eta_val) ? eta_val : NaN)
        end
        if length(md_dense) >= 3
            push!(interp_mdot_rough, _roughness_metric(pr_dense[1:length(md_dense)], md_dense))
        end
        if length(et_dense) >= 3
            push!(interp_eta_rough, _roughness_metric(pr_dense[1:length(et_dense)], et_dense))
        end
    end

    println("Turbine Map Pipeline Diagnostic")
    println("  meanline_path=$(parsed["meanline_path"]) group=$(meanline_group)")
    println("  map_path=$(parsed["map_path"]) group=$(map_group)")
    println("  direct valid samples=$(length(sample_rows))")
    _stats_line("direct roughness mdot(PR)", direct_mdot_rough)
    _stats_line("direct roughness eta(PR)", direct_eta_rough)
    _stats_line("table roughness mdot(PR)", table_mdot_rough)
    _stats_line("table roughness eta(PR)", table_eta_rough)
    _stats_line("interp roughness mdot(PR)", interp_mdot_rough)
    _stats_line("interp roughness eta(PR)", interp_eta_rough)
    _stats_line("map-vs-direct |mdot_corr error|", mdot_err)
    _stats_line("map-vs-direct |eta error|", eta_err)

    if !isempty(direct_eta_rough_by_speed)
        worst_direct = sort(direct_eta_rough_by_speed; by=r -> -r.rough_eta)
        println("  worst direct eta roughness speed lines:")
        for r in Iterators.take(worst_direct, min(5, length(worst_direct)))
            println("    speed_idx=$(r.speed_idx), m_tip=$(r.m_tip), rough_eta=$(r.rough_eta)")
        end
    end
    if !isempty(table_eta_rough_by_speed)
        worst_table = sort(table_eta_rough_by_speed; by=r -> -r.rough_eta)
        println("  worst table eta roughness speed lines:")
        for r in Iterators.take(worst_table, min(5, length(worst_table)))
            println("    speed_idx=$(r.speed_idx), omega_corr=$(r.omega_corr), rough_eta=$(r.rough_eta)")
        end
    end

    d_med = median(_finite_values(direct_eta_rough))
    t_med = median(_finite_values(table_eta_rough))
    i_med = median(_finite_values(interp_eta_rough))
    if isfinite(d_med) && isfinite(t_med) && isfinite(i_med)
        if d_med > 10 * t_med && d_med > 10 * i_med
            println("  diagnosis_hint: direct streamtube samples are noisier than table/interp -> likely meanline/streamtube source")
        elseif t_med > 5 * d_med && t_med > 2 * i_med
            println("  diagnosis_hint: table is noisier than direct/interp -> likely tabulation/root selection source")
        elseif i_med > 5 * t_med && i_med > 5 * d_med
            println("  diagnosis_hint: interpolation is noisier than direct/table -> likely interpolation artifact")
        else
            println("  diagnosis_hint: mixed contributors; inspect per-sample CSV and map generation settings")
        end
    end

    if !isnothing(csv_path)
        path = String(csv_path)
        open(path, "w") do io
            println(io, "speed_idx,m_tip,phi,omega_corr,PR_turb,mdot_corr_direct,eta_direct,map_valid,mdot_corr_map,eta_map,mdot_err,eta_err")
            for r in sample_rows
                println(
                    io,
                    "$(r.speed_idx),$(r.m_tip),$(r.phi),$(r.omega_corr),$(r.PR_turb),$(r.mdot_corr_direct),$(r.eta_direct),$(r.map_valid),$(r.mdot_corr_map),$(r.eta_map),$(r.mdot_err),$(r.eta_err)",
                )
            end
        end
        println("  wrote sample comparison CSV: $(path)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
