# Meanline Model Design

## Goal

Define one turbomachine-level meanline tabulation flow that uses the generic
axial streamtube solver and produces one canonical
`Turbomachine.TabulatedPerformanceMap`.

The intent is to stop generating separate compressor-shaped and turbine-shaped
map formulations at the meanline level. Instead, the meanline model should
sample the same underlying axial machine physics and emit one common map type
that can represent both compression and expansion behavior.

## First Principles

The axial machine streamtube solver already gives us the physically meaningful
machine response for a fixed non-dimensional operating state:
- inlet flow state
- shaft-speed proxy
- axial machine geometry
- blade-row aerodynamic models

From that solve we recover the machine-level outputs we care about:
- `PR = pt_out / pt_in`
- `eta`
- validity / stall / choke information

Those are exactly the outputs needed to build the canonical
`TabulatedPerformanceMap` used by the new operating-point solver.

So the meanline tabulator should be a thin layer that:
1. chooses a tabulation grid
2. converts each dimensional operating point into streamtube inputs
3. samples `streamtube_solve`
4. stores the resulting `PR`, `eta`, and validity region into one shared map

## Canonical Map Relation

The new meanline model should target the same relation for both compressor-like
and turbine-like machines:
- inputs: `omega`, `mdot`
- outputs: `PR`, `eta`

This is the correct target because:
- it matches the new `solve_operating_point(...)` API
- it is symmetric for compressors and turbines
- compressor behavior is the region where `PR > 1`
- turbine behavior is the region where `PR < 1`

The meanline generator should therefore not produce a turbine-specific inverse
map such as `(omega, PR) -> mdot`. That inversion belongs in the operating-point
solver when needed, not in the canonical map construction.

## Canonical Coordinates

The shared `TabulatedPerformanceMap` should be tabulated on corrected
coordinates internally:
- `omega_corr`
- `mdot_corr`

with tables:
- `PR_table[omega_corr_i, mdot_corr_j]`
- `eta_table[omega_corr_i, mdot_corr_j]`

and admissible-flow limits:
- `mdot_corr_min[omega_corr_i]`
- `mdot_corr_max[omega_corr_i]`

This is the correct internal representation because it preserves the
industry-standard corrected-coordinate scaling already used by the legacy
compressor and turbine map implementations.

Callers should continue to work in physical `omega` and physical `mdot`; the
conversion to corrected coordinates belongs inside `TabulatedPerformanceMap`
evaluation and domain queries.

## Proposed Public API

```julia
tabulate_axial_machine_model(
    model::AxialMachine.AxialMachineModel,
    eos::Fluids.AbstractEOS;
    Tt_in_ref::Real,
    pt_in_ref::Real,
    omega_corr_grid,
    mdot_corr_grid,
    Vtheta_inlet::Real=0.0,
    streamtube_radii=AxialMachine.meanline_radii(model),
    prefer_root::Symbol=:low,
    want_diagnostics::Bool=true,
) -> TabulatedPerformanceMap
```

### Inputs
- `model`: axial machine row stack and geometry
- `eos`: thermodynamic model used for dimensionalization
- `Tt_in_ref`, `pt_in_ref`: inlet reference state used for corrected-coordinate
  conversion and dimensional streamtube sampling
- `omega_corr_grid`: corrected-speed grid to tabulate over
- `mdot_corr_grid`: corrected-flow grid to tabulate over
- `Vtheta_inlet`: inlet tangential velocity for the sampled streamtube
- `streamtube_radii`: per-row streamtube radii to sample; defaults to meanline
  radii from the machine geometry
- `prefer_root`: low/high root choice for row-local streamtube branch selection
- `want_diagnostics`: whether to retain per-sample diagnostic payloads

### Outputs
A `TabulatedPerformanceMap` containing:
- `omega_corr` grid
- `mdot_corr` grid
- `PR` table
- `eta` table
- valid flow bounds at each speed
- optional diagnostic metadata if the map type or attached payload supports it

## Proposed Internal Structure

The shared meanline implementation should live in:
- `src/physics/turbomachine/meanline_model.jl`

with the following structure.

### 1. `tabulate_axial_machine_model(...)`
Public entrypoint.

Responsibilities:
- validate inputs
- resolve defaults
- call the helper phases below
- construct the final `TabulatedPerformanceMap`

### 2. `_resolve_meanline_tabulation_grids(...)`
Construct or normalize:
- `omega_corr_grid`
- `mdot_corr_grid`
- `streamtube_radii`

Responsibilities:
- check monotonicity and positivity where required
- ensure streamtube radii are valid for every row
- centralize grid/default handling

### 3. `_sample_axial_machine_streamtube(...)`
Evaluate the streamtube solver at one dimensional operating point.

Inputs:
- `model`
- `eos`
- `Tt_in_ref`, `pt_in_ref`
- `omega`
- `mdot`
- `Vtheta_inlet`
- `streamtube_radii`

Responsibilities:
- convert corrected-coordinate samples to dimensional streamtube inputs
  - `omega_corr -> omega`
  - `mdot_corr -> mdot`
- call `streamtube_solve(...)`
- return:
  - `PR`
  - `eta`
  - validity flags
  - detailed solver diagnostics if requested

### 4. `_compute_feasible_flow_limits(...)`
Determine the valid `mdot_corr` interval at each `omega_corr`.

Responsibilities:
- identify where streamtube samples are valid
- derive `mdot_corr_min` / `mdot_corr_max` per speed line
- optionally refine boundaries with 1D root finding if the coarse grid is not
  sufficient

This is the meanline analogue of the old compressor surge/choke boundary search,
except the implementation should be phrased in generic validity terms rather
than compressor-specific language.

### 5. `_sample_meanline_tables(...)`
Populate the performance tables over the resolved `(omega_corr, mdot_corr)` grid.

Responsibilities:
- call `_sample_axial_machine_streamtube(...)` for each grid point
- store `PR_table`, `eta_table`
- store validity mask / diagnostic payloads if wanted

### 6. `_build_tabulated_performance_map(...)`
Construct the final `TabulatedPerformanceMap` from the sampled data.

Responsibilities:
- package grids and sampled tables into the canonical map object
- attach admissible `mdot` limits
- optionally attach metadata describing reference state and generation options

## Streamtube-Solver Coupling

The generic streamtube solver already operates on the common axial-machine
physics. The new meanline generator should treat it as the source of truth.

The only coupling logic the meanline layer should own is the conversion between
its dimensional tabulation coordinates and the streamtube solver's internal
non-dimensional inputs.

That coupling should remain explicit and localized in one helper rather than
spread throughout separate compressor and turbine tabulation code.

## Compressor And Turbine Interpretation

The new meanline generator should not branch on `machine_kind` to produce two
mathematically different map structures.

Instead:
- a compressor-like axial machine will naturally populate the region where
  `PR > 1`
- a turbine-like axial machine will naturally populate the region where
  `PR < 1`
- a machine capable of both behaviors should still be representable by the same
  map type

This keeps the meanline model aligned with the new common operating-point
solver and avoids baking old map conventions into the generator.

## What To Reuse From Existing Code

Reuse from the current codebase:
- the generic axial streamtube solver in
  `src/physics/turbomachine/axial_machine/streamtube_solver.jl`
- geometry helpers such as `meanline_radii(model)`
- any generic validation or sampling patterns from the existing compressor and
  turbine meanline implementations
- interpolation / table-construction utilities already used by
  `TabulatedPerformanceMap`

Do not carry forward into the canonical path:
- compressor-specific non-dimensional map generation logic
- turbine-specific inverse tabulation logic `(omega, PR) -> mdot`
- wrapper-specific output object construction for compressor or turbine map
  types

## Diagnostics

The shared meanline tabulator should retain enough diagnostics to explain why
particular regions of the sampled domain are invalid.

Useful per-sample diagnostics include:
- raw `streamtube_solve(...)` result
- validity flag
- stall / choke flags
- per-row invalidity markers if available
- the corrected and physical inputs used for that sample

These diagnostics should not be required for the final runtime map, but they are
useful while validating the new common meanline implementation.

## Immediate Refactor Direction

The next implementation should proceed in this order:

1. Add a dedicated design-compatible entrypoint in
   `src/physics/turbomachine/meanline_model.jl` that returns
   `TabulatedPerformanceMap`.
2. Factor the dimensional-to-streamtube conversion into one helper.
3. Implement generic feasible-flow-limit detection in `(omega, mdot)` space.
4. Sample `PR` and `eta` tables directly from the streamtube solver.
5. Construct a canonical `TabulatedPerformanceMap`.
6. Leave the legacy compressor and turbine meanline generators in place only as
   temporary adapters until callers are migrated.

## Practical Consequence For Current Code

The existing turbomachine-level meanline wrapper still delegates to separate
compressor and turbine tabulators and returns wrapper-specific map objects.
That is not the right long-term architecture.

The new meanline generator should instead:
- own the canonical map construction at the turbomachine level
- use the generic axial-machine streamtube solver directly
- emit `TabulatedPerformanceMap`
- let any compressor/turbine-specific formats exist only as compatibility or
  export layers
