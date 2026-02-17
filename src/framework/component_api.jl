abstract type AbstractComponent end

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
