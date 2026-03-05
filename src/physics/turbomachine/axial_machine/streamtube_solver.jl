function _station_rhs(
    gamma::Real,
    pi::Real,
    A::Real,
    tau::Real,
    nu_x::Real,
    nu_theta::Real,
)
    tau > 0 || return NaN
    term = 1 - ((gamma - 1) / (2 * tau)) * (nu_x^2 + nu_theta^2)
    term > 0 || return NaN
    return pi * A * (nu_x / tau) * term^(1 / (gamma - 1))
end

function _bisect_root(
    f::Function,
    a::Float64,
    b::Float64;
    tol::Float64=1e-10,
    max_iters::Int=80,
)
    fa = f(a)
    fb = f(b)
    isfinite(fa) && isfinite(fb) || return NaN
    fa * fb <= 0 || return NaN
    lo = a
    hi = b
    flo = fa
    mid = 0.5 * (lo + hi)
    for _ in 1:max_iters
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        isfinite(fmid) || return NaN
        if abs(fmid) <= tol || abs(hi - lo) <= tol
            return mid
        end
        if flo * fmid <= 0
            hi = mid
        else
            lo = mid
            flo = fmid
        end
    end
    return mid
end

function _find_scalar_roots(
    f::Function,
    x_lo::Float64,
    x_hi::Float64;
    n_scan::Int=201,
    tol::Float64=1e-10,
)
    x_hi > x_lo || return Float64[]
    grid = collect(range(x_lo, x_hi, length=n_scan))
    vals = [f(x) for x in grid]
    roots = Float64[]
    for i in eachindex(grid)
        fi = vals[i]
        if isfinite(fi) && abs(fi) <= tol
            push!(roots, grid[i])
        end
    end
    for i in 1:(length(grid) - 1)
        f1 = vals[i]
        f2 = vals[i + 1]
        if !(isfinite(f1) && isfinite(f2))
            continue
        end
        if f1 * f2 < 0
            r = _bisect_root(f, grid[i], grid[i + 1]; tol=tol)
            isfinite(r) && push!(roots, r)
        end
    end
    isempty(roots) && return roots
    sort!(roots)
    out = Float64[roots[1]]
    for i in 2:length(roots)
        abs(roots[i] - out[end]) > 1e-8 && push!(out, roots[i])
    end
    return out
end

function _solve_station_nu_x(
    gamma::Float64,
    mu::Float64,
    pi::Float64,
    A::Float64,
    tau::Float64,
    nu_theta::Float64;
    prefer::Symbol=:low,
)
    tau > 0 || return (converged=false, nu_x=NaN)
    term = (2 * tau / (gamma - 1)) - nu_theta^2
    term > 0 || return (converged=false, nu_x=NaN)
    x_hi = sqrt(term) * (1 - 1e-8)
    x_lo = 1e-10
    f = nu_x -> mu - _station_rhs(gamma, pi, A, tau, nu_x, nu_theta)
    roots = _find_scalar_roots(f, x_lo, x_hi)
    isempty(roots) && return (converged=false, nu_x=NaN)
    nu_x = prefer == :high ? last(roots) : first(roots)
    return (converged=true, nu_x=nu_x)
end

function nu_u_for_row(
    model::AxialMachineModel,
    row::AxialRow,
    m_tip::Float64,
)
    idx_ref = first_rotor_index(model)
    r_tip_ref = model.rows[idx_ref].r_tip
    if row.kind == :rotor
        return Float64(row.omega_sign) * m_tip * row.r_mean / r_tip_ref
    else
        return 0.0
    end
end

"""
    streamtube_solve(model, m_tip, phi_in; prefer_root=:low)

Run the axial row-marching solve in non-dimensional coordinates.

Inputs:
- `model::AxialMachineModel`: row/station geometry and aero closures.
- `m_tip`: non-dimensional reference rotor-tip speed
  (`m_tip = omega * r_tip_ref / a0_in`).
- `phi_in`: inlet flow coefficient at station 1
  (`phi_in = Vx_1 / U_m1`).
- `prefer_root`: root-selection mode (`:low` or `:high`) when multiple
  continuity roots exist at a station.

Returns a named tuple with:
- `PR`: total-pressure ratio (`Pt_out / Pt_in`), or `NaN` if invalid.
- `eta`: adiabatic efficiency proxy from `(PR, tau_out)`, or `NaN` if invalid.
- `stall::Bool`: true if any row reports stall/invalid aero margin.
- `choke::Bool`: true if any station continuity solve fails.
- `valid::Bool`: true when marching completed and outputs are finite/physical.
- `mu`: conserved corrected massflow-like invariant used by the station solve.
- `tau`, `pi`, `nu_theta`, `nu_x`: per-station state histories.
- `stall_row`, `choke_row`, `valid_row`: per-row diagnostic flags.
"""
function streamtube_solve(
    model::AxialMachineModel,
    m_tip::Real,
    phi_in::Real;
    prefer_root::Symbol=:low,
)
    m_tip_f = Float64(m_tip)
    phi_in_f = Float64(phi_in)
    n_rows = length(model.rows)
    n_stations = n_rows + 1

    tau = fill(NaN, n_stations)
    pi = fill(NaN, n_stations)
    nu_theta = fill(NaN, n_stations)
    nu_x = fill(NaN, n_stations)
    stall_row = falses(n_rows)
    choke_row = falses(n_rows)
    valid_row = trues(n_rows)

    tau[1] = 1.0
    pi[1] = 1.0
    nu_theta[1] = model.nu_theta_station1
    nu_u1 = nu_u_for_row(model, model.rows[first_rotor_index(model)], m_tip_f)
    abs(nu_u1) > 0 || return (
        PR=NaN,
        eta=NaN,
        stall=true,
        choke=true,
        valid=false,
        mu=NaN,
        tau=tau,
        pi=pi,
        nu_theta=nu_theta,
        nu_x=nu_x,
        stall_row=stall_row,
        choke_row=choke_row,
        valid_row=valid_row,
    )
    nu_x[1] = phi_in_f * abs(nu_u1)
    mu = _station_rhs(model.gamma, 1.0, model.A_station[1], 1.0, nu_x[1], nu_theta[1])
    isfinite(mu) || return (
        PR=NaN,
        eta=NaN,
        stall=true,
        choke=true,
        valid=false,
        mu=NaN,
        tau=tau,
        pi=pi,
        nu_theta=nu_theta,
        nu_x=nu_x,
        stall_row=stall_row,
        choke_row=choke_row,
        valid_row=valid_row,
    )

    for k in 1:n_rows
        row = model.rows[k]
        station_in = k
        station_out = k + 1

        inlet = _solve_station_nu_x(
            model.gamma,
            mu,
            pi[station_in],
            model.A_station[station_in],
            tau[station_in],
            nu_theta[station_in];
            prefer=prefer_root,
        )
        if !inlet.converged
            choke_row[k] = true
            valid_row[k] = false
            break
        end
        nu_x[station_in] = inlet.nu_x

        nu_u = nu_u_for_row(model, row, m_tip_f)
        aero_in = RowAeroInput(
            nu_x[station_in],
            nu_theta[station_in],
            nu_u,
            tau[station_in],
            pi[station_in],
            mu,
        )
        aero_out = row_aero(row.aero, row, aero_in)
        stall_row[k] = (aero_out.stall_margin <= 0) || !aero_out.valid
        if !aero_out.valid
            valid_row[k] = false
            break
        end

        gamma_ratio = model.gamma / (model.gamma - 1)
        f_exit = function (nu_x_out)
            if row.kind == :rotor
                nu_theta_out = nu_u + nu_x_out * aero_out.k_theta_exit
                tau_out = tau[station_in] + (model.gamma - 1) * nu_u * (nu_theta_out - nu_theta[station_in])
                tau_out > 0 || return NaN
                pi_out = pi[station_in] * (tau_out / tau[station_in])^gamma_ratio *
                         exp(-gamma_ratio * aero_out.delta_s_hat)
                rhs = _station_rhs(
                    model.gamma,
                    pi_out,
                    model.A_station[station_out],
                    tau_out,
                    nu_x_out,
                    nu_theta_out,
                )
                return mu - rhs
            else
                nu_theta_out = nu_x_out * aero_out.k_theta_exit
                tau_out = tau[station_in]
                pi_out = pi[station_in] * exp(-gamma_ratio * aero_out.delta_s_hat)
                rhs = _station_rhs(
                    model.gamma,
                    pi_out,
                    model.A_station[station_out],
                    tau_out,
                    nu_x_out,
                    nu_theta_out,
                )
                return mu - rhs
            end
        end

        roots = _find_scalar_roots(f_exit, 1e-10, 2.5)
        if isempty(roots)
            choke_row[k] = true
            valid_row[k] = false
            break
        end
        nu_x_out = prefer_root == :high ? last(roots) : first(roots)

        if row.kind == :rotor
            nu_theta[station_out] = nu_u + nu_x_out * aero_out.k_theta_exit
            tau[station_out] = tau[station_in] + (model.gamma - 1) * nu_u *
                               (nu_theta[station_out] - nu_theta[station_in])
            tau[station_out] > 0 || begin
                valid_row[k] = false
                break
            end
            pi[station_out] = pi[station_in] *
                              (tau[station_out] / tau[station_in])^gamma_ratio *
                              exp(-gamma_ratio * aero_out.delta_s_hat)
        else
            nu_theta[station_out] = nu_x_out * aero_out.k_theta_exit
            tau[station_out] = tau[station_in]
            pi[station_out] = pi[station_in] * exp(-gamma_ratio * aero_out.delta_s_hat)
        end
        nu_x[station_out] = nu_x_out
    end

    model_valid = all(valid_row)
    model_choke = any(choke_row)
    model_stall = any(stall_row)
    if !model_valid
        return (
            PR=NaN,
            eta=NaN,
            stall=model_stall,
            choke=model_choke,
            valid=false,
            mu=mu,
            tau=tau,
            pi=pi,
            nu_theta=nu_theta,
            nu_x=nu_x,
            stall_row=stall_row,
            choke_row=choke_row,
            valid_row=valid_row,
        )
    end

    tau_out = tau[end]
    pi_out = pi[end]
    eta_tt = if (tau_out > 0) && (pi_out > 0) && (abs(tau_out - 1) > 1e-12)
        eta = (pi_out^((model.gamma - 1) / model.gamma) - 1) / (tau_out - 1)
        isfinite(eta) ? eta : NaN
    else
        NaN
    end
    return (
        PR=pi_out,
        eta=eta_tt,
        stall=model_stall,
        choke=model_choke,
        valid=model_valid && isfinite(pi_out) && (pi_out > 0) && isfinite(eta_tt),
        mu=mu,
        tau=tau,
        pi=pi,
        nu_theta=nu_theta,
        nu_x=nu_x,
        stall_row=stall_row,
        choke_row=choke_row,
        valid_row=valid_row,
    )
end
