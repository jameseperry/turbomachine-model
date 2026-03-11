# Operating-Point Solver Design

## Goal

Define a single operating-point formulation that works for both compressors and turbines in a network-based turbomachine solver.

The intent is to separate:
- the physical operating-point problem we want to solve
- the representation used by the canonical turbomachine performance map
- the numerical strategy used to solve or sweep that problem

For the new solver, the canonical map type is `TabulatedPerformanceMap`. The same operating-point solver should be used for both compressor and turbine behavior when backed by this map type.

## First Principles

An operating-point solver should answer:

Given fixed inlet conditions, outlet backpressure, and shaft speed, what operating state does the machine settle to?

For both compressors and turbines, the physically meaningful operating-point outputs are:
- `mdot`
- `ht_out`
- `tau`
- `PR = pt_out / pt_in`
- `eta`
- `power = tau * omega`
- all physically admissible branches
- validity / convergence / diagnostics

These are the machine-level quantities that matter to a network solver.

## Canonical Single-Point Problem

### Fixed variables

For both compressors and turbines, the natural fixed inputs are:
- `pt_in`
- `Tt_in` or equivalently `ht_in`
- `pt_out`
- `omega`

### Solved variables

The machine should solve for:
- `mdot`
- `ht_out`
- `tau`

### Derived variables

From the solved state we then report:
- `PR = pt_out / pt_in`
- `eta`
- `power = tau * omega`

## Minimal Unknown Structure

Even though the full physical operating point includes `mdot`, `ht_out`, and `tau`, a map-backed machine usually has a smaller irreducible unknown set.

For `TabulatedPerformanceMap`:
- once `mdot` is known at fixed `(omega, pt_in, Tt_in)`
- the map gives `PR` and `eta`
- `ht_out` follows from the isentropic efficiency relation
- `tau` follows from shaft power balance

So, for many map formulations, the true irreducible unknown is just:
- `mdot`

That means the canonical point-solver interface can still expose `(mdot, ht_out, tau)` as outputs while using a lower-dimensional internal solve.

## Compressor And Turbine Symmetry

The recommended machine-facing operating-point interface is the same for both compressors and turbines:
- fixed: `pt_in`, `ht_in`, `pt_out`, `omega`
- solve: all physically admissible `mdot` branches
- reconstruct: `ht_out`, `tau`
- report: `PR`, `eta`, `power`, diagnostics

With `TabulatedPerformanceMap` as the canonical map type, compressor-like and turbine-like behavior are not different solver APIs. They are different regions of the same map/operating domain, typically distinguished by whether `PR` is above or below 1.0.

This is the cleanest interface for a network-based solver because it reflects the coupling variables between components.

## Sweep Definition

A sweep is not a different physical problem. It is repeated solution of the same operating-point problem over a family of fixed inputs.

A sweep should therefore produce rows of operating points containing:
- sweep coordinates
- solved `mdot`
- solved `ht_out`
- solved `tau`
- derived `PR`
- derived `eta`
- derived `power`
- validity / convergence / diagnostics
- branch identifier, when multiple solutions exist

## Recommended Sweep Variables

For a common compressor/turbine framework, the preferred sweep variables are:
- `omega`
- `pt_out`, or equivalently `PR`

This keeps the sweep symmetric between compressors and turbines and keeps `mdot` as the main solved unknown.

### Recommended 1D sweep
- vary `omega`
- fix `pt_in`, `ht_in`, `pt_out`
- solve `mdot`

### Recommended 2D sweep
- vary `omega`
- vary `pt_out` or `PR`
- fix `pt_in`, `ht_in`
- solve `mdot` at each sweep point

## What Should Not Be Canonical

The following are useful implementation details for map generation or diagnostics, but should not be treated as the canonical public operating-point interface:
- compressor-specific formulations that sweep `omega` and solve for `mdot` only at a target `PR`
- turbine-specific formulations that invert a map in a different direction than the compressor
- formulations based directly on internal meanline variables such as `phi` or `nu_x`

Those may remain useful internally, but they should sit below a common machine-level operating-point API.

## Proposed Canonical API

### Single point

Inputs:
- `map`
- `eos`
- `pt_in`
- `ht_in`
- `pt_out`
- `omega`
- `mdot_guess`

Outputs:
- all physically admissible branches of:
  - `mdot`
  - `ht_out`
  - `tau`
  - `PR`
  - `eta`
  - `power`
- `converged`
- `retcode`
- `diagnostics`

### Sweep

Inputs:
- `map`
- `eos`
- fixed inlet state
- a 1D or 2D sweep over `omega` and `pt_out` or `PR`
- branch-tracking options

Outputs:
- a table of operating points with one row per sweep point or branch

## Numerical Architecture

The point-solver and sweep-solver layers should be separated.

### Point solver responsibilities
- solve one operating point
- return all physically admissible branches and diagnostics
- filter only obviously non-physical candidates

### Sweep solver responsibilities
- choose sweep coordinates
- seed each point solve from nearby prior solutions
- track branches when multiple solutions exist
- handle backoff or feasibility policies if requested
- collect results into arrays or row tables

## Immediate Refactor Direction

A new operating-point solver implementation should be designed around this structure:

1. Define one canonical machine-level point-solver API for both compressors and turbines.
2. Make the point solver solve for the same machine outputs for both machine types.
3. Allow the internal solve to be map-specific.
   - compressor may use scalar `mdot` root-finding
   - turbine may use direct map evaluation or its own scalar solve
4. Build sweep code as orchestration on top of that point solver.
5. Keep branch tracking and continuation in the sweep layer, but make the point solver return every admissible branch.

## Practical Consequence For Current Code

The current code only shared a shallow result-finalization helper between compressor and turbine point solves. That was not a meaningful common operating-point solver architecture.

The next implementation should instead aim for:
- one canonical operating-point problem definition
- one canonical point-solver interface
- one canonical `TabulatedPerformanceMap`-based solver path

That gives a genuinely merged design rather than a thin wrapper around separate implementations.

## API Sketch

This section sketches a concrete public/internal API for the next operating-point solver iteration.

### Public Point Solver

```julia
solve_operating_point(
    map::TabulatedPerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    omega::Real,
    want_diagnostics::Bool=true,
    options::NamedTuple=NamedTuple(),
)
```

#### Inputs
- `map`: canonical turbomachine `TabulatedPerformanceMap`
- `eos`: thermodynamic equation of state
- `pt_in`: inlet stagnation pressure
- `ht_in`: inlet stagnation enthalpy
- `pt_out`: outlet stagnation pressure
- `omega`: shaft speed
- `want_diagnostics`: whether to compute and return extra diagnostics payloads
- `options`: machine-specific or algorithm-specific options bundle

#### Physical admissibility
The point solver should return all branches that are mathematically valid and also
physically admissible according to a minimal filter.

Examples of candidates that should be filtered out:
- negative stagnation pressures
- negative absolute temperatures implied by the thermodynamic reconstruction
- negative stagnation enthalpies if the EOS/model does not permit them
- non-finite values (`NaN`, `Inf`)

Examples of candidates that should not be filtered out merely because they are unusual:
- low-efficiency branches
- windmilling / turbine-like compressor branches
- reverse-torque branches, if the thermodynamics and map remain internally consistent
- multiple branches at the same boundary conditions

```julia
(
    converged=Bool,
    retcode=Symbol,
    candidates=Vector{NamedTuple},
    diagnostics=NamedTuple,
)
```

Where each candidate has at least:

```julia
(
    branch_coordinate=Float64,
    mdot=Float64,
    ht_out=Float64,
    tau=Float64,
    pressure_ratio=Float64,
    efficiency=Float64,
    power=Float64,
    residuals=NTuple{3,Float64},
    physically_admissible=Bool,
    machine_payload=NamedTuple,
)
```

The result may also include optional helper fields such as:

```julia
(
    n_candidates=Int,
    n_admissible_candidates=Int,
)
```

### Public Sweep Solver

```julia
solve_operating_sweep(
    map::TabulatedPerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    omega_grid::AbstractVector,
    pt_out_grid::Union{Nothing,AbstractVector}=nothing,
    pr_grid::Union{Nothing,AbstractVector}=nothing,
    branch_selection::Symbol=:nearest,
    track_branches::Bool=false,
    want_diagnostics::Bool=true,
    options::NamedTuple=NamedTuple(),
)
```

#### Sweep rules
- exactly one of `pt_out_grid` or `pr_grid` should be provided
- if `pr_grid` is provided, each `pt_out` is computed as `pt_in * pr`
- each sweep point calls the point solver
- prior successful branch coordinates are used as continuation seeds
- if `track_branches=true`, all candidates are retained and branch matching is performed across sweep points
- branch selection for single-branch output is applied in the sweep layer, not in the point solver

#### Return value for a single-branch sweep

```julia
(
    mode=:single,
    omega_grid=Vector{Float64},
    pt_out_grid=Vector{Float64},
    mdot=Matrix{Float64},
    ht_out=Matrix{Float64},
    tau=Matrix{Float64},
    PR=Matrix{Float64},
    eta=Matrix{Float64},
    power=Matrix{Float64},
    converged=BitMatrix,
    diagnostics=Matrix{NamedTuple},
)
```

Where `branch_selection` could be:
- `:nearest`
- `:low`
- `:high`
- another explicit branch-selection policy defined by the sweep implementation

#### Return value for an all-branches sweep

```julia
(
    mode=:all,
    rows=Vector{NamedTuple},
    tracking=NamedTuple,
    diagnostics=Matrix{NamedTuple},
)
```

Where each row contains:

```julia
(
    i_omega=Int,
    i_pt_out=Int,
    omega=Float64,
    pt_out=Float64,
    branch_id=Int,
    mdot=Float64,
    ht_out=Float64,
    tau=Float64,
    PR=Float64,
    eta=Float64,
    power=Float64,
    converged=Bool,
    machine_payload=NamedTuple,
)
```

### Internal Adapter API

The public APIs above should be backed by machine-specific adapters that produce candidate operating states.

```julia
_operating_point_candidates(
    map::TabulatedPerformanceMap,
    eos;
    pt_in::Float64,
    ht_in::Float64,
    pt_out::Float64,
    omega::Float64,
    want_diagnostics::Bool,
    options::NamedTuple,
)
```

#### Internal adapter return value

```julia
(
    candidates=Vector{NamedTuple},
    diagnostics=NamedTuple,
)
```

### Canonical Internal Strategy

The canonical `TabulatedPerformanceMap` adapter is expected to:
- compute `target_pr = pt_out / pt_in`
- find one or more `mdot` roots satisfying the map pressure-ratio relation
- reconstruct `ht_out` from efficiency and the isentropic enthalpy change
- reconstruct `tau` from power balance
- evaluate residuals for diagnostics
- emit one candidate per physically admissible branch

Whether a branch behaves like a compressor or a turbine is determined by the solved state, not by dispatching to separate operating-point solver implementations.

### Compatibility Wrappers

During migration, keep these thin wrappers:

```julia
solve_compressor_operating_point(args...; kwargs...) =
    solve_operating_point(args...; kwargs...)

solve_turbine_operating_point(args...; kwargs...) =
    solve_operating_point(args...; kwargs...)
```

and similarly for sweeps if needed.

This allows the codebase to migrate to a common API without immediately breaking existing scripts.

### Optional Helper API

If callers still want a single branch, provide a separate helper instead of making
the point solver itself prune branches:

```julia
select_operating_point_branch(
    point_result;
    policy::Symbol=:nearest,
    reference_coordinate::Union{Nothing,Real}=nothing,
)
```

This keeps the core solver honest: solve first, then select.

## Building Sweep On Top Of Point Solve

The sweep layer should be implemented as orchestration on top of the point solver.

At each sweep point:
1. define the fixed machine boundary conditions for that point
2. call `solve_operating_point(...)`
3. receive all physically admissible branches
4. either:
   - retain all branches and perform branch tracking, or
   - select one branch using a sweep-layer policy
5. store the resulting operating-point data in sweep outputs

This keeps the point solver focused on physics and local admissible branches, while the sweep layer handles continuation and presentation.

### Recommended Point-Solver Extension

To support sweeping well, the point solver should accept optional continuation hints:

```julia
solve_operating_point(
    map::TabulatedPerformanceMap,
    eos::Fluids.AbstractEOS;
    pt_in::Real,
    ht_in::Real,
    pt_out::Real,
    omega::Real,
    continuation_hints::AbstractVector{<:Real}=Float64[],
    want_diagnostics::Bool=true,
    options::NamedTuple=NamedTuple(),
)
```

#### Meaning of `continuation_hints`
- A vector of branch coordinates observed at nearby sweep points
- Used to help the point solver recover the same physical branches consistently
- Especially useful for root-bracketing/root-polishing methods that benefit from prior root locations

For example:
- the canonical map adapter can use them as prior `mdot` roots when scanning for `PR(omega, mdot) = target`

This is a clean extension because it communicates numerical continuation information without polluting the physics contract.

### Candidate Contract For Sweeps

Each point-solver candidate should explicitly include a branch coordinate:

```julia
(
    branch_coordinate=Float64,
    mdot=Float64,
    ht_out=Float64,
    tau=Float64,
    pressure_ratio=Float64,
    efficiency=Float64,
    power=Float64,
    residuals=NTuple{3,Float64},
    physically_admissible=Bool,
    machine_payload=NamedTuple,
)
```

For the current map formulations, `branch_coordinate` will usually just be `mdot`.

Making that explicit avoids hard-coding sweep logic to a particular field name.

## Sweep Algorithm Sketch

### Single-branch sweep

A single-branch sweep selects one branch at each point after the point solve returns all admissible branches.

```julia
prior_branch_coordinates = Float64[]

for point in sweep_grid
    result = solve_operating_point(
        map,
        eos;
        pt_in=point.pt_in,
        ht_in=point.ht_in,
        pt_out=point.pt_out,
        omega=point.omega,
        continuation_hints=prior_branch_coordinates,
        want_diagnostics=want_diagnostics,
        options=options,
    )

    selected = select_operating_point_branch(
        result;
        policy=branch_selection,
        reference_coordinate=previous_selected_branch_coordinate,
    )

    store selected output
    prior_branch_coordinates = [c.branch_coordinate for c in result.candidates]
end
```

The sweep layer is responsible for choosing the branch-selection policy, not the point solver.

### All-branches sweep

An all-branches sweep stores every admissible branch returned by the point solver, then optionally matches them across neighboring sweep points.

```julia
prior_branch_coordinates = Float64[]
all_candidates = Vector{Vector{NamedTuple}}()

for point in sweep_grid
    result = solve_operating_point(
        map,
        eos;
        pt_in=point.pt_in,
        ht_in=point.ht_in,
        pt_out=point.pt_out,
        omega=point.omega,
        continuation_hints=prior_branch_coordinates,
        want_diagnostics=want_diagnostics,
        options=options,
    )

    push!(all_candidates, result.candidates)
    prior_branch_coordinates = [c.branch_coordinate for c in result.candidates]
end

tracking = track_branches_over_grid(all_candidates)
```

This is where branch tracking belongs: in the sweep layer, after the point solver has honestly returned all admissible local branches.

## Sweep-Layer Policies

The sweep layer may apply policies that should not live inside the point solver itself.

Examples:
- `branch_selection = :low`
- `branch_selection = :high`
- `branch_selection = :nearest`
- backoff of requested `pt_out` or `PR` if a point is infeasible
- interpolation of missing isolated sweep points for plotting only
- branch-matching cost functions for continuation/tracking

These are sweep-orchestration concerns, not point-physics concerns.

## Suggested Helper Functions

To support sweep-on-top-of-point cleanly, the implementation should likely expose these helpers:

```julia
select_operating_point_branch(
    point_result;
    policy::Symbol=:nearest,
    reference_coordinate::Union{Nothing,Real}=nothing,
)
```

```julia
branch_coordinates(point_result) -> Vector{Float64}
```

```julia
track_operating_point_branches(
    point_results::AbstractVector;
    distance::Function,
    max_match_cost::Float64=Inf,
)
```

These keep the point solver small while making the sweep layer easy to implement.

## Recommended Division Of Responsibility

### Point solver
- solve one operating point
- return all physically admissible local branches
- accept continuation hints
- report diagnostics for that local solve

### Sweep solver
- define sweep coordinates
- call point solver repeatedly
- propagate continuation hints across neighboring points
- select a branch, or keep all branches
- track branches over the sweep grid
- apply backoff or presentation-oriented policies if requested
- collect results into structured arrays or row tables

## Consequence For Refactoring

This design suggests the following refactor direction:

1. Implement a branch-complete `solve_operating_point(...)`
2. Add `continuation_hints` to its API
3. Add an explicit `branch_coordinate` field to each candidate
4. Rewrite `solve_operating_sweep(...)` to call `solve_operating_point(...)` at each sweep point
5. Move branch selection and branch tracking fully into the sweep layer

That would give a real layered architecture:
- local physics solve at the point level
- global continuation/tracking logic at the sweep level
