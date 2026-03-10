"""
Generic operating-point point-solve orchestration.

This utility captures the common solver skeleton shared by compressor and turbine
point solves:
- fixed unknown vector shape `(mdot, ht_out, tau)`,
- nonlinear solve on a 3-equation residual,
- optional scaled residuals for numerical conditioning,
- optional diagnostics payload returned by the residual callback.

Physics stays outside this file via callbacks.
"""

using NonlinearSolve

Base.@kwdef struct OperatingPointSolveOptions
    # If true, the nonlinear loop uses `residual_eval(...).scaled`.
    # If false, the loop uses `residual_eval(...).unscaled`.
    scaled_residuals::Bool = true

    # Diagnostics policy:
    # - :never      -> never request diagnostics payload
    # - :on_failure -> request diagnostics only after a failed solve
    # - :always     -> request diagnostics after every solve
    diagnostics_mode::Symbol = :on_failure

    abstol::Float64 = 1e-10
    reltol::Float64 = 1e-8
    maxiters::Int = 100
end

"""
    solve_operating_point_generic(
        context;
        u0,
        residual_eval!,
        project_solution,
        build_scales=nothing,
        options=OperatingPointSolveOptions(),
    )

Generic 3-unknown operating-point solver.

Expected unknown ordering is fixed to:
- `u[1] = mdot`
- `u[2] = ht_out`
- `u[3] = tau`

This matches the current compressor/turbine point-solve formulations.

Callback contract:

1. `build_scales(context, u0_vector) -> scales_payload` (optional)
   Called once before solve. Return any payload needed by residual scaling.
   If omitted, `scales_payload === nothing`.

2. `residual_eval!(R, u_vector, context, scales_payload; scaled::Bool, want_diagnostics::Bool)`
   Required in-place callback that writes residuals into `R`.

   Requirements:
   - `R` is length-3 and must be overwritten in-place.
   - `scaled=true` means write scaled residuals.
   - `scaled=false` means write unscaled residuals.
   - if `want_diagnostics=false`, return value is ignored (recommend `nothing`).
   - if `want_diagnostics=true`, return an arbitrary diagnostics `NamedTuple`.

   Important usage pattern:
   - During nonlinear iterations, this callback is invoked with
     `want_diagnostics=false` to keep the inner loop cheap.
   - After solve, it is invoked once more for final unscaled residuals and
     optional diagnostics.

3. `project_solution(u_vector, context) -> payload`
   Required callback for post-solve projection into model-specific outputs
   (for example map values, derived power terms, operating mode labels, etc.).

Returns a named tuple:
- `mdot, ht_out, tau` : solved unknowns
- `u`                 : raw solved vector `[mdot, ht_out, tau]`
- `projected`         : payload from `project_solution`
- `residuals_unscaled`: final unscaled residual vector
- `diagnostics`       : diagnostics payload or `nothing`
- `scales`            : scales payload or `nothing`
- `retcode, converged, solution`
"""
function solve_operating_point_generic(
    context;
    u0::AbstractVector{<:Real},
    residual_eval!::Function,
    project_solution::Function,
    build_scales::Union{Nothing,Function}=nothing,
    options::OperatingPointSolveOptions=OperatingPointSolveOptions(),
)
    options.abstol > 0 || error("abstol must be > 0")
    options.reltol > 0 || error("reltol must be > 0")
    options.maxiters >= 1 || error("maxiters must be >= 1")
    options.diagnostics_mode in (:never, :on_failure, :always) ||
        error("diagnostics_mode must be one of: :never, :on_failure, :always")

    u0_vec = _to_u0_vector(u0)
    scales = isnothing(build_scales) ? nothing : build_scales(context, u0_vec)

    function f!(R, u, _)
        # Inner loop: do not request diagnostics to keep iteration cost down.
        residual_eval!(
            R,
            u,
            context,
            scales;
            scaled=options.scaled_residuals,
            want_diagnostics=false,
        )
        _validate_written_residuals!(R, "residual_eval! output")
        return nothing
    end

    prob = NonlinearProblem(f!, u0_vec)
    sol = solve(
        prob,
        NewtonRaphson(; autodiff=AutoForwardDiff());
        abstol=options.abstol,
        reltol=options.reltol,
        maxiters=options.maxiters,
    )

    converged = (string(sol.retcode) == "Success")
    u_sol = Float64.(sol.u)

    # Always compute final unscaled residuals for reporting/validation.
    residuals_unscaled_vec = zeros(3)
    residual_eval!(
        residuals_unscaled_vec,
        u_sol,
        context,
        scales;
        scaled=false,
        want_diagnostics=false,
    )
    _validate_written_residuals!(residuals_unscaled_vec, "final unscaled residuals")
    residuals_unscaled = (
        residuals_unscaled_vec[1],
        residuals_unscaled_vec[2],
        residuals_unscaled_vec[3],
    )

    want_diag = (
        options.diagnostics_mode == :always ||
        (options.diagnostics_mode == :on_failure && !converged)
    )
    diagnostics = if want_diag
        R_diag = zeros(3)
        diag = residual_eval!(
            R_diag,
            u_sol,
            context,
            scales;
            scaled=false,
            want_diagnostics=true,
        )
        _validate_written_residuals!(R_diag, "diagnostics residual output")
        diag isa NamedTuple || error(
            "residual_eval!(...; want_diagnostics=true) must return a NamedTuple, got $(typeof(diag))",
        )
        diag
    else
        nothing
    end

    projected = project_solution(u_sol, context)

    return (
        mdot=u_sol[1],
        ht_out=u_sol[2],
        tau=u_sol[3],
        u=u_sol,
        projected=projected,
        residuals_unscaled=residuals_unscaled,
        diagnostics=diagnostics,
        scales=scales,
        retcode=sol.retcode,
        converged=converged,
        solution=sol,
    )
end

function _to_u0_vector(u0::AbstractVector{<:Real})
    length(u0) == 3 || error("u0 must have length 3 in [mdot, ht_out, tau] order")
    return Float64[u0[1], u0[2], u0[3]]
end

function _validate_written_residuals!(R, label::AbstractString)
    length(R) == 3 || error("$(label) must have length 3")
    for i in 1:3
        v = _primal_value(R[i])
        isfinite(v) || error("$(label) contains non-finite value at index $(i): $(v)")
    end
    return nothing
end
