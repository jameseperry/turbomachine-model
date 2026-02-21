# TurboMachineModel Component-First Schema Sketch (v3)

This document defines a component-first internal schema for engine modeling in Julia.

Core idea: you instantiate concrete component types (`CompressorSection`, `Combustor`, `TurbineSection`, etc.) directly. Parameters are constructor arguments for each component type, not fields on a generic `ComponentInstance` record.

Current implementation status: the live code currently exposes only a minimal structural network API in `src/network.jl` (`Network`, `EndpointRef`, `add_component!`, `connect!`). The rest of this document is forward-looking design.

## 1. Design Principles

1. Concrete component types carry behavior.
2. Multiple dispatch provides the component interface.
3. Network graph stores component objects and connectivity.
4. Ports are acausal by default; connections define direct variable coupling constraints.
5. Assembly builds a global ODE/DAE from component + connection equations.

## 2. Core Schema

### 2.1 Component Hierarchy

```julia
abstract type AbstractComponent end

struct CompressorSection <: AbstractComponent
    map_id::Symbol
    eta_guess::Float64
    init::NamedTuple
end

struct Combustor <: AbstractComponent
    dp_frac::Float64
    fuel_lhv::Float64
    init::NamedTuple
end

struct TurbineSection <: AbstractComponent
    map_id::Symbol
    eta_guess::Float64
    init::NamedTuple
end

struct ShaftInertia <: AbstractComponent
    J::Float64
    damping::Float64
    init::NamedTuple
end
```

### 2.2 Ports and Endpoints

```julia
struct EndpointRef
    component::Symbol
    port::Symbol
end
```

Port definitions are provided by component methods, for example:

- `port_specs(c::CompressorSection)`
- `port_specs(c::Combustor)`
- `port_specs(c::TurbineSection)`
- `port_specs(c::ShaftInertia)`

Each port spec declares at least:

1. `domain` (`:fluid`, `:shaft`, `:thermal`, `:signal`)
2. `role` (`:interface`, `:node`, `:control`)
3. `variables` (for example `[:p,:T,:h,:mdot,:composition]`)

### 2.3 Network Graph

```julia
struct Network
    components::Dict{Symbol,AbstractComponent}
    # each: (from, to, vars)
    # vars can be :all or an explicit vector of variable symbols
    connections::Vector{NamedTuple{
        (:from, :to, :vars),
        Tuple{EndpointRef,EndpointRef,Union{Symbol,Vector{Symbol}}}
    }}
    # each: (id, target, bc_type, value_spec, priority)
    boundaries::Vector{NamedTuple{
        (:id, :target, :bc_type, :value_spec, :priority),
        Tuple{Symbol,EndpointRef,Symbol,Any,Int}
    }}
end
```

### 2.4 Assembly-Time Controls

Assembly and solver controls are passed into functions instead of being stored in `Network`.

Examples:

- `assemble(model; formulation=:dae, index_reduction=:auto, diagnostics=:basic)`
- `solve_system(sys, tspan; initialization=:steady_guess, reltol=1e-6, abstol=1e-8)`

## 3. Component API (Draft)

The interface is defined by generic functions and dispatch, not a keyword-based trait system.

### 3.1 Structural / Model-Building Functions

These functions are used to construct and validate a legal model graph.

Required:

1. `port_specs(c)::Dict{Symbol,PortSpec}`

Recommended:

1. `parameter_symbols(c)::Vector{Symbol}`  # for diagnostics/config tooling

### 3.2 Equation / Integration Support Functions

These functions define unknowns and equations used by assembly and integration.

Required:

1. `state_symbols(c, node_id)::Vector{Symbol}`  # differential unknowns owned by component
2. `algebraic_symbols(c, node_id)::Vector{Symbol}`  # algebraic unknowns owned by component (empty allowed)
3. `residuals!(res, c, ctx)::Nothing`  # writes component residual equations into preallocated `res`
4. `initial_conditions(c, node_id)::Dict{Symbol,Any}`  # initial values for component-owned states

### 3.3 Optional Performance / Runtime Hooks

1. `jacobian_sparsity(c, ctx)`  # sparse pattern contribution
2. `jacobian_values!(Jv, c, ctx)`  # analytic Jacobian contribution
3. `events(c, ctx)`  # event surfaces / switching support
4. `output_symbols(c, node_id)` and `outputs!(y, c, ctx)`  # logging/telemetry

### 3.4 Minimal `ctx` Contract for `residuals!`

`ctx` should provide at least:

1. Access to current values of component-owned unknowns.
2. Access to connected port variable values.
3. Access to simulation time and constant model parameters.

## 4. Port Typing and Compatibility (Draft)

Port typing is explicit so connection validation can be deterministic.
`Unitful.jl` is used for unit metadata and compatibility checks.
Connections are ideal interfaces in this model, so transported properties are equality-coupled.
Connections are strictly 2-port only (`connect!` enforces max one connection per endpoint).

```julia
using Unitful

struct PortSpec
    domain::Symbol                 # :fluid, :shaft, :thermal, :signal
    # entries like (var=:p, unit=u"Pa", coupling=:equality)
    # coupling is :equality (X_from = X_to) or :conservation (X_from + X_to = 0)
    vars::Vector{NamedTuple{(:var, :unit, :coupling), Tuple{Symbol,Unitful.Units,Symbol}}}
end
```

Each component declares typed ports through:

- `port_specs(c)::Dict{Symbol,PortSpec}`

Example intent:

- Fluid port might declare:
  - `vars = [(var=:p, unit=u"Pa", coupling=:equality), (var=:T, unit=u"K", coupling=:equality), (var=:h, unit=u"J/kg", coupling=:equality), (var=:mdot, unit=u"kg/s", coupling=:conservation), (var=:composition, unit=u"1", coupling=:equality)]`
- Shaft port might declare:
  - `vars = [(var=:omega, unit=u"rad/s", coupling=:equality), (var=:tau, unit=u"N*m", coupling=:conservation)]`

### 4.1 Compatibility Rules

For each connection `(from, to, vars)`:

1. Resolve `from_spec` and `to_spec` from `port_specs`.
2. Require `from_spec.domain == to_spec.domain`.
3. Require the connection links exactly two ports (one `from`, one `to`).
4. Determine coupled vars:
   - if `vars == :all`, use set intersection policy only if exact sets match for now.
   - otherwise `vars` is explicit and every variable must exist on both ports.
5. For each coupled variable, require matching `coupling` type.
6. For `:equality`, add `X_from = X_to`.
7. For `:conservation`, add `X_from + X_to = 0` using outward-positive convention on both ports.
8. Require unit convertibility for each coupled variable (`uconvert`-compatible).
9. (Optional, strict mode) require exact same units in addition to convertibility.

### 4.2 Coupling Equation Intent

Given a validated variable:

1. `:equality` -> add equality constraint (`x_from = x_to`)
2. `:conservation` -> add anti-equality constraint (`x_from + x_to = 0`)
3. Through variables use outward-positive sign convention on all ports.

### 4.3 Validation Pseudocode

```julia
function validate_connection(model, conn)
    from_spec = port_specs(model.components[conn.from[1]])[conn.from[2]]
    to_spec   = port_specs(model.components[conn.to[1]])[conn.to[2]]
    from_vars = Dict(v.var => (unit=v.unit, coupling=v.coupling) for v in from_spec.vars)
    to_vars   = Dict(v.var => (unit=v.unit, coupling=v.coupling) for v in to_spec.vars)

    from_spec.domain == to_spec.domain || error("domain mismatch")

    used_vars = conn.vars == :all ? collect(keys(from_vars)) : conn.vars
    conn.vars == :all && Set(keys(from_vars)) == Set(keys(to_vars)) || error("vars mismatch for :all")

    for v in used_vars
        v in keys(from_vars) || error("missing var on source: $v")
        v in keys(to_vars)   || error("missing var on target: $v")
        from_vars[v].coupling == to_vars[v].coupling || error("coupling mismatch on $v")
        (from_vars[v].coupling in (:equality, :conservation)) || error("invalid coupling on $v")
        uconvert(1.0 * from_vars[v].unit, 1.0 * to_vars[v].unit)  # throws if incompatible
    end
end
```

## 5. Variable Conventions (MVP)

- Fluid ports: `[:p, :T, :h, :mdot, :composition]`
- Shaft ports: `[:omega, :tau]`
- Thermal ports: `[:T, :qdot]`

## 6. Port Coupling Examples (Draft)

These are concrete examples of port variable sets for turbomachinery models.

### 6.1 Fluid Through Port (compressor/turbine/duct interface)

- `(var=:p, unit=u"Pa", coupling=:equality)`
- `(var=:T, unit=u"K", coupling=:equality)`
- `(var=:h, unit=u"J/kg", coupling=:equality)`
- `(var=:mdot, unit=u"kg/s", coupling=:conservation)`
- `(var=:composition, unit=u"1", coupling=:equality)`

### 6.2 Shaft Port

- `(var=:omega, unit=u"rad/s", coupling=:equality)`
- `(var=:tau, unit=u"N*m", coupling=:conservation)`

### 6.3 Thermal Port

- `(var=:T_wall, unit=u"K", coupling=:equality)`
- `(var=:qdot, unit=u"W", coupling=:conservation)`

### 6.4 Fuel Injection Port

- `(var=:p_fuel, unit=u"Pa", coupling=:equality)`
- `(var=:T_fuel, unit=u"K", coupling=:equality)`
- `(var=:mdot_fuel, unit=u"kg/s", coupling=:conservation)`
- `(var=:fuel_type, unit=u"1", coupling=:equality)`  # could also be composition scalar/vector

### 6.5 Bleed / Bypass Fluid Port

- `(var=:p, unit=u"Pa", coupling=:equality)`
- `(var=:T, unit=u"K", coupling=:equality)`
- `(var=:h, unit=u"J/kg", coupling=:equality)`
- `(var=:mdot, unit=u"kg/s", coupling=:conservation)`
- `(var=:composition, unit=u"1", coupling=:equality)`

### 6.6 Generator Load Port (mechanical side)

- `(var=:omega, unit=u"rad/s", coupling=:equality)`
- `(var=:tau_load, unit=u"N*m", coupling=:conservation)`

### 6.7 Command Variables

Command interfaces are future work and not part of the current minimal `Network` API.

## 7. Example Usage Sketch

This is illustrative data-shape pseudocode, not finalized implementation code.

```julia
cmp = CompressorSection(
    map_id=:cmp_map_a,
    eta_guess=0.86,
    init=(plenum_mass=0.8,),
)
cmb = Combustor(
    dp_frac=0.04,
    fuel_lhv=43e6,
    init=(internal_energy=1.2e6,),
)
trb = TurbineSection(
    map_id=:trb_map_a,
    eta_guess=0.90,
    init=(plenum_mass=0.6,),
)
spool = ShaftInertia(
    J=0.35,
    damping=0.01,
    init=(omega=800.0,),
)

network = Network(
    components = Dict(
        :cmp   => cmp,
        :cmb   => cmb,
        :trb   => trb,
        :spool => spool,
    ),

    connections = [
        (from=EndpointRef(:cmp,:outlet), to=EndpointRef(:cmb,:inlet),   vars=:all),
        (from=EndpointRef(:cmb,:outlet), to=EndpointRef(:trb,:inlet),   vars=:all),
        (from=EndpointRef(:cmp,:shaft),  to=EndpointRef(:spool,:left),  vars=[:omega, :tau]),
        (from=EndpointRef(:trb,:shaft),  to=EndpointRef(:spool,:right), vars=[:omega, :tau]),
    ],

    boundaries = [
        (id=:b1, target=EndpointRef(:cmp,:inlet),   bc_type=:fixed_total_pressure,    value_spec=101325u"Pa", priority=10),
        (id=:b2, target=EndpointRef(:cmp,:inlet),   bc_type=:fixed_total_temperature, value_spec=288.15u"K",  priority=10),
        (id=:b3, target=EndpointRef(:trb,:outlet),  bc_type=:fixed_static_pressure,   value_spec=101325u"Pa", priority=10),
        (id=:b4, target=(:spool,:right), bc_type=:applied_torque_load,     value_spec=0.0u"N*m",   priority=5),
    ),
)
```

## 8. Assembly Flow

1. Validate component constructor invariants.
2. Validate referenced ports against `port_specs(c)`.
3. Validate connection compatibility.
4. Build global variable index from component states and connection variables.
5. Assemble local component equations via `equations!`.
6. Add connection equations by coupling type (`:equality -> X_from = X_to`, `:conservation -> X_from + X_to = 0`).
7. Gather initial conditions from each component and form final ODE/DAE problem data (assembly/solver options are function arguments).

## 9. Notes

- This schema removes generic `type_id` / `equation_provider_id` indirection.
- Connection behavior is intentionally minimal: strict 2-port ideal interfaces with per-variable coupling type.
- Through variables use outward-positive sign convention.
- Assembly controls are intentionally not stored on `Network`.
- Initial conditions are component-owned (`init`) rather than centralized.
- Multi-port behavior (tees/manifolds/mixing/losses) should be modeled as explicit components, not connections.
- Extension path: add a new component type and implement required methods.
- Typed API, builder API, and DSL can all target this same internal model representation.
