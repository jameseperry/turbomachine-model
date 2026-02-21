module Component

import ..Port: port_shape

abstract type AbstractComponent end

const PortMapping = NamedTuple{(:port,:label),Tuple{Symbol,Symbol}}

"""
Canonical component variable identity.
`id` is unique within the component scope.
"""
struct ComponentVariable
    id::Symbol
    unit::Any
    mappings::Vector{PortMapping}
end

function ComponentVariable(;
    id::Symbol,
    unit,
    mappings::Vector{PortMapping}=PortMapping[],
)
    seen = Set{Tuple{Symbol,Symbol}}()
    for m in mappings
        key = (m.port, m.label)
        key in seen && error("duplicate mapping ($(m.port), $(m.label)) for variable :$id")
        push!(seen, key)
    end
    return ComponentVariable(id, unit, mappings)
end

"""
A concrete component port has a name and references a shared `PortShape`.
"""
struct ComponentPort
    shape_id::Symbol
    name::Symbol
end

function ComponentPort(; shape_id::Symbol, name::Symbol)
    port_shape(shape_id)  # validate shape exists
    return ComponentPort(shape_id, name)
end

"""
Return structural port specifications for a component.
Expected return type: `Vector{ComponentPort}`.
Components must overload this.
"""
function ports(::AbstractComponent)
    error("ports not implemented for this component type")
end

function has_port(c::AbstractComponent, port_name::Symbol)
    for p in ports(c)
        p.name == port_name && return true
    end
    return false
end

"""
Return component variables.
Expected return type: `Vector{ComponentVariable}`.
`sim_type` should be `:steady` or `:transient`.
Components must overload this.
"""
function variables(::AbstractComponent, sim_type::Symbol)
    error("variables not implemented for this component type")
end

export AbstractComponent, ports
export ComponentVariable, ComponentPort
export has_port
export variables

end # module Component
