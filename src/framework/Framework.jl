module Framework

include("component_api.jl")
include("port_presets.jl")
include("types.jl")
include("network_api.jl")
include("topology.jl")
include("validate_network.jl")

export AbstractComponent, port_specs, required_ports, validate
export FLUID_THROUGH_VARS, SHAFT_VARS, THERMAL_VARS
export FLUID_PORT, SHAFT_PORT, THERMAL_PORT
export EndpointRef, Endpoint, VariableSpec, PortSpec
export ConnectionSpec, BoundarySpec, Network, Model
export ValidationIssue, ValidationReport
export add_component!, connect!, add_boundary!, validate_network
export adjacency, connections_for, component_neighbors, domain_connections

"""
Structural modeling framework primitives
(component API, ports, model types, and validation hooks).
"""

end # module Framework
