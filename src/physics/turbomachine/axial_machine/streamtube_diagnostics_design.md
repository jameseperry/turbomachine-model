# Streamtube Diagnostics Design

## Goal
Add enough structured output around `streamtube_solve(...)` to diagnose row-by-row
behavior of an axial machine without changing the core nondimensional solver into
a dimensional one.

This is specifically aimed at distinguishing:

1. a bad demo machine definition from
2. a bad streamtube/meanline model implementation.

## Current State
`streamtube_solve(...)` already computes most of the useful nondimensional state:

- station arrays:
  - `tau`
  - `pi`
  - `nu_theta`
  - `nu_x`
- row validity flags:
  - `stall_row`
  - `choke_row`
  - `valid_row`

The main issue is that the solver does not currently preserve the intermediate
row-level quantities that explain *why* the stage is behaving the way it is.

In particular, `_advance_row!` computes and then discards:

- `nu_u`
- `k_theta_exit`
- `delta_s_hat`
- `stall_margin`
- `incidence`
- `theta_in`
- `theta_out`
- `nu_theta_out - nu_theta_in`
- `tau_out - tau_in`
- `pi_out - pi_in`

Those are the quantities needed to tell whether a rotor is compressing or
extracting work.

## Design Principles

1. Keep the core solver nondimensional.
2. Do not make plotting scripts reconstruct row internals themselves.
3. Separate:
   - raw nondimensional solve output
   - dimensional reconstruction
   - presentation / scripting
4. Preserve solver clarity: diagnostics should be collected, not entangled with
   the numerical marching logic.

## Proposed Changes

### 1. Enrich `streamtube_solve(...)` output

Keep the existing array outputs, but add two structured payloads:

- `station_data::Vector{NamedTuple}`
- `row_data::Vector{NamedTuple}`

### 2. Return row diagnostics from `_advance_row!`

Instead of returning only convergence / stall / choke flags, `_advance_row!`
should also return a diagnostics payload containing the row-local aerodynamic
and thermodynamic transition data.

Suggested return shape:

```julia
(
    converged=true,
    choke=false,
    stall=false,
    diagnostics=(
        row_radius=...,
        nu_u=...,
        nu_x_in=...,
        nu_x_out=...,
        nu_theta_in=...,
        nu_theta_out=...,
        tau_in=...,
        tau_out=...,
        pi_in=...,
        pi_out=...,
        k_theta_exit=...,
        delta_s_hat=...,
        stall_margin=...,
        incidence=...,
        theta_in=...,
        theta_out=...,
    ),
)
```

`incidence`, `theta_in`, and `theta_out` are already available through
`blade_aero(...).diagnostics`.

### 3. Construct `station_data`

At the end of `streamtube_solve(...)`, build a per-station summary from the
arrays that already exist.

Suggested shape:

```julia
station_data[k] = (
    station_index=k,
    area=station_area(model, k),
    nu_x=nu_x[k],
    nu_theta=nu_theta[k],
    tau=tau[k],
    pi=pi[k],
)
```

For row-based radii:

- station 1 can use the first row radius
- station N+1 can use the last row radius
- interior stations can use the average of adjacent row radii

This should be done explicitly in a helper so the convention is visible.

### 4. Construct `row_data`

Collect one row diagnostics entry per solved row:

```julia
row_data[k] = (
    row_index=k,
    speed_ratio_to_ref=row.speed_ratio_to_ref,
    r_hub=row.r_hub,
    r_tip=row.r_tip,
    row_radius=row_radius,
    valid=...,
    stall=...,
    choke=...,
    nu_u=...,
    nu_x_in=...,
    nu_x_out=...,
    nu_theta_in=...,
    nu_theta_out=...,
    delta_nu_theta=nu_theta_out - nu_theta_in,
    tau_in=...,
    tau_out=...,
    delta_tau=tau_out - tau_in,
    pi_in=...,
    pi_out=...,
    delta_pi=pi_out - pi_in,
    k_theta_exit=...,
    delta_s_hat=...,
    stall_margin=...,
    incidence=...,
    theta_in=...,
    theta_out=...,
)
```

This is the minimum useful diagnostic dataset for turbine/compressor sign
analysis.

## Dimensional Reconstruction Layer

Do **not** dimensionalize inside the core solver.

Instead add a separate helper layer, likely in a new file:

- `src/physics/turbomachine/axial_machine/streamtube_diagnostics.jl`

Suggested helper:

```julia
reconstruct_streamtube_physical_states(model, result, eos; pt_in, Tt_in, omega)
```

This helper would convert the nondimensional quantities into dimensional values:

- `Vx`
- `Vtheta`
- `V`
- `Tt`
- `pt`
- `ht`
- optionally static:
  - `T`
  - `p`
  - `rho`
  - `Mach`

### Why separate this?

Because the solver is correctly formulated in nondimensional variables. The
dimensional view is a *postprocessing concern*, not a solver concern.

## Higher-Level Diagnostic API

Add a thin wrapper above the raw solver:

```julia
diagnose_streamtube_operating_point(
    model,
    eos;
    m_tip,
    phi_in,
    nu_theta_inlet=0.0,
    Tt_in,
    pt_in,
    streamtube_radii=meanline_radii(model),
)
```

This should:

1. call `streamtube_solve_with_phi(...)`
2. collect nondimensional `station_data` and `row_data`
3. optionally call the dimensional reconstruction helper
4. return a structured result for scripts

Suggested top-level return shape:

```julia
(
    summary=(PR=..., eta=..., valid=..., stall=..., choke=...),
    station_data=...,
    row_data=...,
    physical_station_data=...,
    physical_row_data=...,
)
```

## Minimum Useful First Pass

The fastest path to useful diagnostics is to add only nondimensional outputs:

### Per row
- `nu_u`
- `nu_theta_in`
- `nu_theta_out`
- `delta_nu_theta`
- `tau_in`
- `tau_out`
- `delta_tau`
- `pi_in`
- `pi_out`
- `delta_pi`
- `incidence`
- `theta_in`
- `theta_out`
- `stall_margin`

### Per station
- `nu_x`
- `nu_theta`
- `tau`
- `pi`
- `area`

### Geometry/context values worth including
- per station:
  - `area`
- per row:
  - `row_annulus_area`
  - `area_in`
  - `area_out`
  - `theta_metal_in`
  - `theta_metal_out`

That is sufficient to answer:

- Is the rotor extracting work or adding it?
- Is the rotor increasing or decreasing absolute swirl?
- Is pressure rising or falling through the machine?
- Are odd maps caused by the demo geometry / aero setup or by the solver?

## Recommended Implementation Order

1. Modify `_advance_row!` to return diagnostics.
2. Extend `streamtube_solve(...)` to assemble `row_data` and `station_data`.
3. Add a thin `diagnose_streamtube_operating_point(...)` helper.
4. Only after that, add dimensional reconstruction helpers.
5. Then build a script on top of the diagnostic API.

## What Not To Do

- Do not redimensionalize the core solver.
- Do not put all reconstruction logic into scripts.
- Do not make the scripts reverse-engineer row internals from the final arrays.
- Do not mix formatted presentation with solver internals.

## Main Payoff

With these changes, one direct diagnostic run will show whether a supposed
"turbine" stage is actually:

- adding swirl and raising `tau` like a compressor, or
- removing swirl and dropping `tau` like a turbine.

That is the most direct way to separate a bad demo definition from a bad
streamtube model.
