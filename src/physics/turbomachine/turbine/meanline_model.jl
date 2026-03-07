"""
Turbine meanline model wrappers built on shared axial-machine AxialMachine.
"""

import ..AxialMachine
using ....Utility: bracket_bisect_roots

const TurbineMeanlineModel = AxialMachine.AxialMachineModel

function _turbine_eta_from_compressor_eta(eta_comp::Float64)
    if eta_comp < -1e-8
        return clamp(-1.0 / eta_comp, 0.0, 1.0)
    end
    return NaN
end

function _turbine_state_from_phi(
    model::TurbineMeanlineModel,
    m_tip::Float64,
    phi_in::Float64,
    rho0_in_ref::Float64,
    inlet_area::Float64,
    mean_radius_inlet::Float64,
    a0_in_ref::Float64,
    Tt_in_ref::Float64,
    Pt_in_ref::Float64,
    Tt_ref::Float64,
    Pt_ref::Float64,
)
    vals = AxialMachine.streamtube_solve_with_phi(model, m_tip, phi_in)
    vals.valid || return (valid=false, PR_turb=NaN, eta=NaN, mdot=NaN, mdot_corr=NaN)
    isfinite(vals.PR) || return (valid=false, PR_turb=NaN, eta=NaN, mdot=NaN, mdot_corr=NaN)
    vals.PR > 0 || return (valid=false, PR_turb=NaN, eta=NaN, mdot=NaN, mdot_corr=NaN)
    vals.PR < 1 || return (valid=false, PR_turb=NaN, eta=NaN, mdot=NaN, mdot_corr=NaN)
    eta_t = _turbine_eta_from_compressor_eta(Float64(vals.eta))
    isfinite(eta_t) || return (valid=false, PR_turb=NaN, eta=NaN, mdot=NaN, mdot_corr=NaN)

    omega = m_tip * a0_in_ref / model.r_tip_ref
    mdot = phi_in * rho0_in_ref * inlet_area * abs(omega * mean_radius_inlet)
    mdot_corr = corrected_flow(mdot, Tt_in_ref, Pt_in_ref, Tt_ref, Pt_ref)
    return (
        valid=true,
        PR_turb=(1.0 / vals.PR),
        eta=eta_t,
        mdot=mdot,
        mdot_corr=mdot_corr,
    )
end

function _resolve_turbine_tabulation_grids(
    model::TurbineMeanlineModel,
    n_speed::Int,
    n_pr::Int,
    omega_corr_grid::Union{Nothing,Vector{<:Real}},
    pr_turb_grid::Union{Nothing,Vector{<:Real}},
    Tt_in_ref::Float64,
    Tt_ref::Float64,
)
    n_speed >= 2 || error("n_speed must be >= 2")
    n_pr >= 2 || error("n_pr must be >= 2")
    m_lo, m_hi = model.m_tip_bounds
    a0_in_ref = sqrt(model.gamma * model.gas_constant * Tt_in_ref)
    omega_lo = m_lo * a0_in_ref / model.r_tip_ref
    omega_hi = m_hi * a0_in_ref / model.r_tip_ref
    omega_corr_lo = corrected_speed(omega_lo, Tt_in_ref, Tt_ref)
    omega_corr_hi = corrected_speed(omega_hi, Tt_in_ref, Tt_ref)
    omega_corr = isnothing(omega_corr_grid) ?
        collect(range(omega_corr_lo, omega_corr_hi, length=n_speed)) :
        Float64.(omega_corr_grid)
    length(omega_corr) >= 2 || error("omega_corr_grid must have at least 2 points")
    issorted(omega_corr) || error("omega_corr_grid must be sorted ascending")

    pr_grid = isnothing(pr_turb_grid) ? nothing : Float64.(pr_turb_grid)
    if !isnothing(pr_grid)
        length(pr_grid) >= 2 || error("pr_turb_grid must have at least 2 points")
        issorted(pr_grid) || error("pr_turb_grid must be sorted ascending")
        all(pr_grid .> 1.0) || error("pr_turb_grid must be > 1.0 for turbine tabulation")
    end

    return (omega_corr_grid=omega_corr, pr_turb_grid=pr_grid, a0_in_ref=a0_in_ref)
end

function _sorted_pr_samples(samples::Vector{NamedTuple})
    sorted = sort(samples; by=s -> s.PR_turb)
    pr = Float64[]
    mdot = Float64[]
    eta = Float64[]
    i = 1
    n = length(sorted)
    while i <= n
        pr_i = Float64(sorted[i].PR_turb)
        mdot_sum = 0.0
        eta_sum = 0.0
        c = 0
        j = i
        while j <= n && Float64(sorted[j].PR_turb) == pr_i
            mdot_sum += Float64(sorted[j].mdot_corr)
            eta_sum += Float64(sorted[j].eta)
            c += 1
            j += 1
        end
        push!(pr, pr_i)
        push!(mdot, mdot_sum / c)
        push!(eta, eta_sum / c)
        i = j
    end
    return (pr=pr, mdot=mdot, eta=eta)
end

function _branch_pr_samples(
    samples::Vector{NamedTuple},
    branch::Symbol;
    merge_rtol::Float64=1e-4,
)
    isempty(samples) && return (pr=Float64[], mdot=Float64[], eta=Float64[])
    branch in (:low, :high) || error("branch must be one of: low|high")
    merge_rtol >= 0 || error("merge_rtol must be >= 0")

    sorted = sort(samples; by=s -> s.PR_turb)
    pr = Float64[]
    mdot = Float64[]
    eta = Float64[]

    i = 1
    n = length(sorted)
    while i <= n
        pr_i = Float64(sorted[i].PR_turb)
        tol = merge_rtol * max(abs(pr_i), 1.0)
        j = i
        candidate_idx = i
        candidate_phi = Float64(sorted[i].phi)
        while j <= n && abs(Float64(sorted[j].PR_turb) - pr_i) <= tol
            phi_j = Float64(sorted[j].phi)
            if branch == :high
                if phi_j > candidate_phi
                    candidate_phi = phi_j
                    candidate_idx = j
                end
            else
                if phi_j < candidate_phi
                    candidate_phi = phi_j
                    candidate_idx = j
                end
            end
            j += 1
        end

        s = sorted[candidate_idx]
        push!(pr, Float64(s.PR_turb))
        push!(mdot, Float64(s.mdot_corr))
        push!(eta, Float64(s.eta))
        i = j
    end

    return (pr=pr, mdot=mdot, eta=eta)
end

function _interp_pr_sample(
    pr_axis::Vector{Float64},
    mdot_axis::Vector{Float64},
    eta_axis::Vector{Float64},
    pr_target::Float64,
)
    n = length(pr_axis)
    n >= 1 || return (mdot_corr=NaN, eta=NaN)
    if pr_target <= pr_axis[1]
        return (mdot_corr=mdot_axis[1], eta=eta_axis[1])
    elseif pr_target >= pr_axis[end]
        return (mdot_corr=mdot_axis[end], eta=eta_axis[end])
    end

    i_hi = searchsortedfirst(pr_axis, pr_target)
    i_hi <= 1 && return (mdot_corr=mdot_axis[1], eta=eta_axis[1])
    i_hi > n && return (mdot_corr=mdot_axis[end], eta=eta_axis[end])
    i_lo = i_hi - 1
    pr_lo = pr_axis[i_lo]
    pr_hi = pr_axis[i_hi]
    if pr_hi == pr_lo
        return (mdot_corr=mdot_axis[i_lo], eta=eta_axis[i_lo])
    end
    t = (pr_target - pr_lo) / (pr_hi - pr_lo)
    mdot_val = (1 - t) * mdot_axis[i_lo] + t * mdot_axis[i_hi]
    eta_val = (1 - t) * eta_axis[i_lo] + t * eta_axis[i_hi]
    return (mdot_corr=mdot_val, eta=eta_val)
end

function _smooth_axis(values::Vector{Float64}, window::Int)
    window <= 1 && return copy(values)
    n = length(values)
    n <= 2 && return copy(values)
    half = window ÷ 2
    out = similar(values)
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        out[i] = sum(@view values[lo:hi]) / (hi - lo + 1)
    end
    return out
end

"""
Tabulate an axial-machine model into a dimensional tabulated turbine map.

The axial solver is sampled in `(m_tip, phi_in)` and converted to turbine map
coordinates `(omega_corr, PR_turb)` where `PR_turb = Pt_in / Pt_out`.
"""
function tabulate_turbine_meanline_model(
    model::TurbineMeanlineModel;
    n_speed::Int=31,
    n_pr::Int=41,
    omega_corr_grid::Union{Nothing,Vector{<:Real}}=nothing,
    pr_turb_grid::Union{Nothing,Vector{<:Real}}=nothing,
    Tt_in_ref::Real=288.15,
    Pt_in_ref::Real=101_325.0,
    Tt_ref::Real=Tt_in_ref,
    Pt_ref::Real=Pt_in_ref,
    interpolation::Symbol=:bilinear,
    boundary_resolution::Int=401,
    branch::Symbol=:high,
    pr_grid_mode::Symbol=:union,
    sampling_strategy::Symbol=:roots,
    pr_merge_rtol::Real=1e-4,
    smooth_window::Int=0,
)
    boundary_resolution >= 21 || error("boundary_resolution must be >= 21")
    interpolation in (:bilinear, :bicubic) ||
        error("interpolation must be :bilinear or :bicubic")
    branch in (:low, :high) || error("branch must be one of: low|high")
    pr_grid_mode in (:union, :overlap) || error("pr_grid_mode must be one of: union|overlap")
    sampling_strategy in (:roots, :scan) || error("sampling_strategy must be one of: roots|scan")
    pr_merge_rtol >= 0 || error("pr_merge_rtol must be >= 0")
    smooth_window >= 0 || error("smooth_window must be >= 0")
    smooth_window <= 1 || isodd(smooth_window) || error("smooth_window must be odd when > 1")
    Tt_in_ref > 0 || error("Tt_in_ref must be > 0")
    Pt_in_ref > 0 || error("Pt_in_ref must be > 0")
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")

    grids = _resolve_turbine_tabulation_grids(
        model,
        n_speed,
        n_pr,
        omega_corr_grid,
        pr_turb_grid,
        Float64(Tt_in_ref),
        Float64(Tt_ref),
    )
    omega_corr = grids.omega_corr_grid
    a0_in_ref = grids.a0_in_ref

    rho0_in_ref = Pt_in_ref / (model.gas_constant * Tt_in_ref)
    inlet_area = AxialMachine.station_area(model, 1)
    mean_radius_inlet = abs(model.speed_ratio_ref) * model.r_flow_ref
    phi_lo, phi_hi = model.phi_in_bounds
    phi_scan = collect(range(phi_lo, phi_hi, length=boundary_resolution))

    valid_speed_idx = Int[]
    m_tip_valid = Float64[]
    omega_corr_valid = Float64[]
    pr_min = Float64[]
    pr_max = Float64[]
    speed_samples = Vector{NamedTuple}()

    for (i, omega_corr_i) in pairs(omega_corr)
        omega_i = omega_corr_i * sqrt(Tt_in_ref / Tt_ref)
        m_tip_i = omega_i * model.r_tip_ref / a0_in_ref
        samples_i = NamedTuple[]
        for phi in phi_scan
            st = _turbine_state_from_phi(
                model,
                m_tip_i,
                phi,
                rho0_in_ref,
                inlet_area,
                mean_radius_inlet,
                a0_in_ref,
                Float64(Tt_in_ref),
                Float64(Pt_in_ref),
                Float64(Tt_ref),
                Float64(Pt_ref),
            )
            st.valid || continue
            push!(samples_i, (phi=phi, PR_turb=st.PR_turb, eta=st.eta, mdot_corr=st.mdot_corr))
        end
        isempty(samples_i) && continue
        push!(valid_speed_idx, i)
        push!(m_tip_valid, m_tip_i)
        push!(omega_corr_valid, omega_corr_i)
        push!(pr_min, minimum(s.PR_turb for s in samples_i))
        push!(pr_max, maximum(s.PR_turb for s in samples_i))
        push!(speed_samples, (m_tip=m_tip_i, omega_corr=omega_corr_i, samples=samples_i))
    end

    length(valid_speed_idx) >= 2 ||
        error("meanline model has fewer than two valid turbine-like speed lines in requested tabulation range")

    pr_grid = if isnothing(grids.pr_turb_grid)
        pr_lo, pr_hi = if pr_grid_mode == :union
            (minimum(pr_min), maximum(pr_max))
        else
            (maximum(pr_min), minimum(pr_max))
        end
        pr_hi > pr_lo || error("invalid turbine PR_turb bounds derived from meanline model for pr_grid_mode=$(pr_grid_mode)")
        collect(range(pr_lo, pr_hi, length=n_pr))
    else
        grids.pr_turb_grid
    end

    mdot_corr_table = Matrix{Float64}(undef, length(m_tip_valid), length(pr_grid))
    eta_table = Matrix{Float64}(undef, length(m_tip_valid), length(pr_grid))

    for (i, speed_data) in pairs(speed_samples)
        m_tip_i = speed_data.m_tip
        samples_i = speed_data.samples
        branch_axis = _branch_pr_samples(
            samples_i,
            branch;
            merge_rtol=Float64(pr_merge_rtol),
        )
        if length(branch_axis.pr) < 2
            branch_axis = _sorted_pr_samples(samples_i)
        end
        prior_phi = Float64[]
        mdot_row = Vector{Float64}(undef, length(pr_grid))
        eta_row = Vector{Float64}(undef, length(pr_grid))
        for (j, pr_target) in pairs(pr_grid)
            sample = if sampling_strategy == :roots
                residual = function (phi)
                    st = _turbine_state_from_phi(
                        model,
                        m_tip_i,
                        phi,
                        rho0_in_ref,
                        inlet_area,
                        mean_radius_inlet,
                        a0_in_ref,
                        Float64(Tt_in_ref),
                        Float64(Pt_in_ref),
                        Float64(Tt_ref),
                        Float64(Pt_ref),
                    )
                    st.valid || return NaN
                    return st.PR_turb - pr_target
                end
                roots = bracket_bisect_roots(
                    residual,
                    (phi_lo, phi_hi);
                    n_scan=boundary_resolution,
                    root_tol=1e-8,
                    prior_roots=prior_phi,
                )

                if !isempty(roots)
                    phi_sel = branch == :high ? last(roots) : first(roots)
                    st = _turbine_state_from_phi(
                        model,
                        m_tip_i,
                        phi_sel,
                        rho0_in_ref,
                        inlet_area,
                        mean_radius_inlet,
                        a0_in_ref,
                        Float64(Tt_in_ref),
                        Float64(Pt_in_ref),
                        Float64(Tt_ref),
                        Float64(Pt_ref),
                    )
                    st.valid || error("root-selected turbine sample invalid at speed=$(speed_data.omega_corr), PR_turb=$(pr_target)")
                    prior_phi = [phi_sel]
                    (mdot_corr=st.mdot_corr, eta=st.eta)
                else
                    # Fallback when root bracketing misses/disconnects: use smooth
                    # 1D interpolation on scanned PR->(mdot_corr, eta) for this speed line.
                    _interp_pr_sample(
                        branch_axis.pr,
                        branch_axis.mdot,
                        branch_axis.eta,
                        Float64(pr_target),
                    )
                end
            else
                # Scan strategy: tabulate directly from branch-selected scan samples.
                _interp_pr_sample(
                    branch_axis.pr,
                    branch_axis.mdot,
                    branch_axis.eta,
                    Float64(pr_target),
                )
            end

            mdot_row[j] = sample.mdot_corr
            eta_row[j] = sample.eta
        end

        if smooth_window > 1
            mdot_row = _smooth_axis(mdot_row, smooth_window)
            eta_row = _smooth_axis(eta_row, smooth_window)
        end
        mdot_corr_table[i, :] = mdot_row
        eta_table[i, :] = eta_row
    end

    return TabulatedTurbinePerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        omega_corr_valid,
        pr_grid,
        mdot_corr_table,
        eta_table;
        interpolation=interpolation,
        pr_turb_min=pr_min,
        pr_turb_max=pr_max,
    )
end

"""
Demo turbine meanline model for development/testing.
"""
demo_turbine_meanline_model() = AxialMachine.demo_axial_turbine_model()
