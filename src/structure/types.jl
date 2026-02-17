const Endpoint = Tuple{Symbol,Symbol}

const VariableSpec = NamedTuple{
    (:var, :unit, :coupling),
    Tuple{Symbol,Any,Symbol},
}

struct PortSpec
    domain::Symbol
    vars::Vector{VariableSpec}
end

const ConnectionSpec = NamedTuple{
    (:from, :to, :vars),
    Tuple{Endpoint,Endpoint,Union{Symbol,Vector{Symbol}}},
}

const BoundarySpec = NamedTuple{
    (:id, :target, :bc_type, :value_spec, :priority),
    Tuple{Symbol,Endpoint,Symbol,Any,Int},
}

struct Model
    components::Dict{Symbol,AbstractComponent}
    connections::Vector{ConnectionSpec}
    boundaries::Vector{BoundarySpec}
end
