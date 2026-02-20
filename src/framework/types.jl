struct EndpointRef
    component::Symbol
    port::Symbol
end

const Endpoint = EndpointRef

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
    Tuple{EndpointRef,EndpointRef,Union{Symbol,Vector{Symbol}}},
}

const BoundarySpec = NamedTuple{
    (:id, :target, :bc_type, :value_spec, :priority),
    Tuple{Symbol,EndpointRef,Symbol,Any,Int},
}

struct Network
    components::Dict{Symbol,AbstractComponent}
    connections::Vector{ConnectionSpec}
    boundaries::Vector{BoundarySpec}
end

const Model = Network

struct ValidationIssue
    level::Symbol
    code::Symbol
    message::String
end

struct ValidationReport
    valid::Bool
    issues::Vector{ValidationIssue}
end
