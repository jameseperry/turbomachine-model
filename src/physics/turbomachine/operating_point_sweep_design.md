# Operating-Point Sweep Script Design

## Goal

Define a single generic operating-point sweep script for turbomachines using:
- `TabulatedPerformanceMap`
- `solve_operating_point(...)`
- `solve_operating_sweep(...)`

The script should work for both compressor-like and turbine-like operating regions.

## Purpose Of The Script

The sweep script is a diagnostic and visualization tool. Its job is to:
1. define a family of fixed operating conditions
2. solve the operating point over that family
3. visualize the solved machine response and feasibility structure

It is not a separate physics model.

## Canonical Sweep Variables

The script should sweep over:
- `omega`
- `PR = pt_out / pt_in`

Equivalent alternative:
- `omega`
- `pt_out`

Recommendation:
- expose both in the CLI
- normalize internally to `pt_out`
- compute `PR` for plotting and reporting

This is symmetric for compressors and turbines and matches the network-facing machine problem.

## Fixed Inputs

The script should require or default:
- `pt_in`
- `Tt_in` or `ht_in`
- `eos`
- input map path

Recommended default inlet state:
- `pt_in = 101325.0`
- `Tt_in = 288.15`
- `eos = ideal air`

## Canonical Solver Call

At each sweep point, the script should call:

```julia
solve_operating_sweep(
    map,
    eos;
    pt_in=...,
    ht_in=...,
    omega_grid=...,
    pt_out_grid=...,
    branch_selection=...,
    track_all_branches=...,
    initial_branch_coordinate=...,
    options=...,
)
```

The script should not reimplement root finding, continuation, or branch tracking itself.

## Output Data Model

### Single-branch mode

Used when the user wants one selected branch per sweep point.

Store matrix-valued fields:
- `omega_grid`
- `pt_out_grid`
- `PR_grid`
- `mdot`
- `ht_out`
- `tau`
- `power`
- `eta`
- `converged`
- `branch_coordinate`
- `branch_index`
- `diagnostics`

### All-branches mode

Used when the user wants every admissible branch.

Store row data:
- `i_omega`
- `i_pt_out`
- `omega`
- `pt_out`
- `PR`
- `branch_id`
- `branch_coordinate`
- `mdot`
- `ht_out`
- `tau`
- `power`
- `eta`
- `converged`
- `machine_payload`

Also retain:
- tracking metadata
- per-point diagnostics

## Plot Set

The same plot set should be used for compressors and turbines.

### Core plots

1. `mdot` over `(omega, PR)`
- shows the flow response of the machine
- likely the most important operating-point plot

2. `eta` over `(omega, PR)`
- shows where the machine is efficient or poor

3. `power` over `(omega, PR)`
- captures compressor power consumption and turbine power extraction in one field
- sign of power is physically meaningful

4. feasibility / branch-count over `(omega, PR)`
- shows where the solver found no branch, one branch, or multiple branches

### Optional secondary plots

5. `tau` over `(omega, PR)`
- often useful if torque coupling matters directly

6. branch-coordinate plot over `(omega, PR)`
- usually `mdot`
- redundant in single-branch mode, useful in all-branches diagnostics

## Compressor And Turbine Defaults

The script should share one implementation but offer different default sweep ranges.

### Compressor defaults
- `PR_min = 1.0`
- `PR_max > 1.0`
- branch selection default: `:high` or `:nearest`

### Turbine defaults
- `PR_min < 1.0`
- `PR_max = 1.0`
- branch selection default: `:nearest`

These are defaults only. The script should not enforce machine kind by physics dispatch.

## Recommended CLI

```text
julia plot_operating_sweep.jl <map-path> \
    --pt-in 101325 \
    --Tt-in 288.15 \
    --omega-min 600 \
    --omega-max 1000 \
    --n-omega 41 \
    --pr-min 1.0 \
    --pr-max 2.5 \
    --n-pr 41 \
    --branch-selection nearest \
    --track-all-branches false \
    --machine-kind compressor \
    --output-prefix out/plot
```

### CLI options
- positional:
  - `map-path`
- inlet state:
  - `--pt-in`
  - `--Tt-in`
- omega sweep:
  - `--omega-min`
  - `--omega-max`
  - `--n-omega`
- pressure loading sweep:
  - either:
    - `--pr-min`, `--pr-max`, `--n-pr`
  - or:
    - `--pt-out-min`, `--pt-out-max`, `--n-pt-out`
- branch policy:
  - `--branch-selection low|high|nearest`
  - `--initial-branch-coordinate <value>`
  - `--track-all-branches`
- presentation defaults:
  - `--machine-kind compressor|turbine`
- outputs:
  - `--output-prefix`
  - `--write-csv`
  - `--write-plots`

## Machine Kind

`machine-kind` should only control defaults and labeling.

It should not change the solver algorithm.

Examples:
- choose default PR range
- choose plot titles
- choose whether to describe positive power as "consumed" or "extracted"

## CSV Output

The script should always be able to write a row-oriented CSV.

### Single-branch CSV schema
- `i_omega`
- `i_pt_out`
- `omega`
- `pt_out`
- `PR`
- `mdot`
- `ht_out`
- `tau`
- `power`
- `eta`
- `converged`
- `branch_coordinate`
- `branch_index`

### All-branches CSV schema
- all of the above, plus:
- `branch_id`
- `local_candidate_index`

## Plotting Recommendations

### Single-branch mode
Use heatmaps or contours for:
- `mdot`
- `eta`
- `power`
- branch count / feasibility

### All-branches mode
Use:
- feasibility / branch-count heatmap
- scatter or line overlays for branch traces
- optional separate branch-id panels if useful

The key point is that all-branches data is not naturally a single smooth scalar field. Do not force it into a contour plot when multiple branches coexist.

## Suggested Implementation Structure

```text
scripts/plot_operating_sweep.jl
```

Internal structure:
1. parse CLI
2. load map
3. build inlet state
4. build `omega_grid` and `pt_out_grid`
5. call `solve_operating_sweep(...)`
6. convert results to row/matrix forms as needed
7. write CSV
8. generate plots

Helper functions to keep in the script or a nearby module:
- `_resolve_sweep_grids(...)`
- `_rows_from_single_branch_sweep(...)`
- `_rows_from_all_branch_sweep(...)`
- `_plot_single_branch_sweep(...)`
- `_plot_all_branch_sweep(...)`
- `_default_pr_range(machine_kind)`

## What To Do With Existing Scripts

The likely end state should be:
- replace `plot_compressor_operating_sweep.jl`
- replace `plot_turbine_operating_sweep.jl`
- keep thin wrappers only if needed for CLI compatibility

Preferred final interface:
- one script: `plot_operating_sweep.jl`
- two wrapper scripts at most, only to preserve current demo commands

## Recommended Refactor Order

1. implement the generic plotting script
2. make it work for compressor defaults
3. make it work for turbine defaults
4. optionally add thin compatibility wrappers:
   - `plot_compressor_operating_sweep.jl`
   - `plot_turbine_operating_sweep.jl`
5. update `demo/Makefile`

## Decision Summary

The correct first-principles sweep script is:
- one generic script
- sweep over `omega` and `PR` or `pt_out`
- use the same solver for compressor and turbine regions
- treat compressor/turbine as labeling/defaults, not different algorithms
- show `mdot`, `eta`, `power`, and feasibility as the primary outputs
