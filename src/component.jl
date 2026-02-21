module Component

import ..Port: port_shape

abstract type AbstractComponent end

"""
Canonical component variable identity.
`id` is unique within the component scope.
"""
struct ComponentVariable
    id::Symbol
    unit::Any
end

"""
A concrete component port is a mapping from a shared `PortShape` to
component-specific variable identities.
"""
struct ComponentPort
    shape_id::Symbol
    variable_ids::Dict{Symbol,Symbol}
end

function ComponentPort(shape_id::Symbol; variable_ids::Union{Nothing,Dict{Symbol,Symbol}}=nothing)
    shape = port_shape(shape_id)
    ids = Dict(v.label => v.label for v in shape.variables)

    if variable_ids !== nothing
        for (label, variable_id) in variable_ids
            haskey(ids, label) ||
                error("unknown label :$label for port shape $shape_id")
            ids[label] = variable_id
        end
    end

    return ComponentPort(shape_id, ids)
end

"""
Return structural port specifications for a component.
Expected return type: `Dict{Symbol,ComponentPort}`.
Components must overload this.
"""
function ports(::AbstractComponent)
    error("ports not implemented for this component type")
end

"""
Return component variables.
Expected return type: `Vector{ComponentVariable}`.
Components must overload this.
"""
function variables(::AbstractComponent)
    error("variables not implemented for this component type")
end

export AbstractComponent, ports
export ComponentVariable, ComponentPort
export variables

end # module Component
