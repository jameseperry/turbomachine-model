Network() = Network(
    Dict{Symbol,AbstractComponent}(),
    ConnectionSpec[],
    BoundarySpec[],
)

function _endpoint_exists(net::Network, ep::EndpointRef)
    comp = get(net.components, ep.component, nothing)
    comp === nothing && return false
    return ep.port in keys(port_specs(comp))
end

function add_component!(net::Network, id::Symbol, component::AbstractComponent)
    haskey(net.components, id) &&
        error("component id already exists: $id")
    validate(component)
    net.components[id] = component
    return net
end

function _endpoint_connection_count(net::Network, ep::EndpointRef)
    count(c -> c.from == ep || c.to == ep, net.connections)
end

function connect!(
    net::Network,
    a::EndpointRef,
    b::EndpointRef;
    vars::Union{Symbol,Vector{Symbol}}=:all,
)
    _endpoint_exists(net, a) || error("unknown endpoint: ($(a.component), $(a.port))")
    _endpoint_exists(net, b) || error("unknown endpoint: ($(b.component), $(b.port))")
    a == b && error("connection endpoints must be distinct")
    _endpoint_connection_count(net, a) == 0 ||
        error("2-port rule violation: endpoint already connected ($(a.component), $(a.port))")
    _endpoint_connection_count(net, b) == 0 ||
        error("2-port rule violation: endpoint already connected ($(b.component), $(b.port))")

    push!(net.connections, (from=a, to=b, vars=vars))
    return net
end

function add_boundary!(
    net::Network,
    id::Symbol,
    target::EndpointRef,
    bc_type::Symbol,
    value_spec;
    priority::Int=10,
)
    _endpoint_exists(net, target) ||
        error("unknown boundary target endpoint: ($(target.component), $(target.port))")
    any(b -> b.id == id, net.boundaries) &&
        error("boundary id already exists: $id")

    push!(
        net.boundaries,
        (
            id=id,
            target=target,
            bc_type=bc_type,
            value_spec=value_spec,
            priority=priority,
        ),
    )
    return net
end
