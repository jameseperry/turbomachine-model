module Components

include("abstract_component.jl")
include("port_presets.jl")

using .AbstractComponentDef
using .PortPresets

const AbstractComponent = AbstractComponentDef.AbstractComponent

"""
Return the structural port specification for a component.
Components should overload this.
"""
function port_specs(::AbstractComponent)
    error("port_specs not implemented for this component type")
end

"""
Return required port names for a component.
Components should overload this.
"""
function required_ports(::AbstractComponent)
    Symbol[]
end

"""
Validate a component's parameterization and invariants.
Components should overload this.
"""
function validate(c::AbstractComponent)
    c
end

include("fluid/Fluid.jl")
include("mechanical/Mechanical.jl")

using .Fluid
using .Mechanical

export AbstractComponent, PortPresets, Fluid, Mechanical, port_specs, required_ports, validate

"""
Namespace for concrete component implementations
(compressor, turbine, combustor, shaft, etc.).
Structural component data types are defined; behavior is deferred.
"""

end # module Components
