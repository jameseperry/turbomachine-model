module Port

struct PortVariable
    label::Symbol
    role::Symbol
    unit::Any
end

struct PortShape
    variables::Vector{PortVariable}
end

function PortVariable(; label::Symbol, role::Symbol, unit)
    role in (:equality, :conservation) ||
        error("invalid role: $role (expected :equality or :conservation)")
    return PortVariable(label, role, unit)
end

function PortShape(; variables::Vector{PortVariable})
    isempty(variables) && error("port shape must define at least one variable")

    labels = Set{Symbol}()
    for v in variables
        v.label in labels && error("duplicate label :$(v.label) in port shape")
        push!(labels, v.label)
    end

    return PortShape(variables)
end

const PORT_SHAPES = Dict{Symbol,PortShape}()

function register_port_shape!(
    id::Symbol,
    shape::PortShape;
    overwrite::Bool=false,
)
    haskey(PORT_SHAPES, id) && !overwrite &&
        error("port shape already registered: $id")
    PORT_SHAPES[id] = shape
    return shape
end

has_port_shape(id::Symbol) = haskey(PORT_SHAPES, id)

function port_shape(id::Symbol)
    shape = get(PORT_SHAPES, id, nothing)
    shape === nothing && error("unknown port shape: $id")
    return shape
end

registered_port_shapes() = collect(keys(PORT_SHAPES))

register_port_shape!(:FluidPort,
    PortShape(;variables=[
            PortVariable(label=:pt, role=:equality, unit=:Pa),
            PortVariable(label=:ht, role=:equality, unit=:J_per_kg),
            PortVariable(label=:composition, role=:equality, unit=:species),
            PortVariable(label=:mdot, role=:conservation, unit=:kg_per_s),
        ],
    )
);

register_port_shape!(:ShaftPort,
    PortShape(;variables=[
            PortVariable(label=:omega, role=:equality, unit=:rad_per_s),
            PortVariable(label=:tau, role=:conservation, unit=:N_m),
        ],
    )
);

register_port_shape!(:ThermalPort,
    PortShape(;variables=[
            PortVariable(label=:temperature, role=:equality, unit=:K),
            PortVariable(label=:heat_flow, role=:conservation, unit=:W),
        ],
    )
)

export PortVariable, PortShape
export register_port_shape!, has_port_shape, port_shape, registered_port_shapes

end # module Port
