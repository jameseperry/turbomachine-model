# Components Interface

This directory contains concrete component implementations.

## Base Type

All components should subtype:

- `AbstractComponent`

Defined in:

- `src/components/abstract_component.jl`

## Required Structural Interface

These methods are required for structural model assembly.

1. `port_specs(c)::Dict{Symbol,Any}`
- Returns port specs by port name.
- Port specs should follow the shape used by `PortPresets` and `Structure.PortSpec`.

2. `required_ports(c)::Vector{Symbol}`
- Returns the list of required ports for this component.
- Each required port must be connected or boundary-constrained by model validation.

3. `validate(c)::AbstractComponent`
- Validates component parameters and invariants.
- Should throw a descriptive error on invalid configuration.
- Returns `c` when valid (convention used in constructors).

## Optional Structural Interface

1. `boundary_targets(c)::Vector{Tuple{Symbol,Symbol}}`
- Allowed boundary targets as `(port_id, variable)` pairs.

2. `parameter_symbols(c)::Vector{Symbol}`
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

1. Keep port definitions centralized where possible using `PortPresets`.
2. Use outward-positive sign convention for through variables.
3. Keep connection semantics in structure/validation code; keep component files focused on component behavior and invariants.
