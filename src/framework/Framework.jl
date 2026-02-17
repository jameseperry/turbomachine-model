module Framework

include("component_api.jl")
include("port_presets.jl")
include("types.jl")

export AbstractComponent, port_specs, required_ports, validate
export FLUID_THROUGH_VARS, SHAFT_VARS, FLUID_PORT, SHAFT_PORT
export Endpoint, VariableSpec, PortSpec, ConnectionSpec, BoundarySpec, Model

"""
Structural modeling framework primitives
(component API, ports, model types, and validation hooks).
"""

end # module Framework
