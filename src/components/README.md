# Components Interface

This directory contains concrete component implementations.

## Base Type

All components should subtype:

- `AbstractComponent`

Defined in:

- `src/component.jl`

## Required Structural Interface

These methods are required for structural model assembly.

1. `ports(c)::Dict{Symbol,Any}`
- Returns port specs by port name.
- Port specs should follow the shape used by component presets and `Component.ComponentPort`.

## Optional Structural Interface

1. `parameter_symbols(c)::Vector{Symbol}`
- Parameter names for diagnostics/config tooling.

## Optional Integration Interface

These are not required for structural-only assembly, but are expected for simulation.

1. `state_symbols(c, node_id)::Vector{Symbol}`
2. `algebraic_symbols(c, node_id)::Vector{Symbol}`
3. `residuals!(res, c, ctx)::Nothing`
4. `initial_conditions(c, node_id)::Dict{Symbol,Any}`

## Optional Performance Hooks

1. `jacobian_sparsity(c, ctx)`
2. `jacobian_values!(Jv, c, ctx)`
3. `events(c, ctx)`
4. `output_symbols(c, node_id)` and `outputs!(y, c, ctx)`

## Conventions

1. Keep port definitions explicit in each component constructor via `ComponentPort`.
2. Use outward-positive sign convention for through variables.
3. Keep network semantics in `src/network.jl`; keep component files focused on component behavior and invariants.
4. Network topology is managed in `src/network.jl` (`Network`, `EndpointRef`, `add_component!`, `connect!`).
