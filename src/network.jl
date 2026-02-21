module Network

import ..Component: AbstractComponent, has_port

struct EndpointRef
    component::Symbol
    port::Symbol
end

struct Connection
    id::Int
    first::EndpointRef
    second::EndpointRef
end

mutable struct Network
    components::Dict{Symbol,AbstractComponent}
    connections_by_id::Dict{Int,Connection}
    connection_by_endpoint::Dict{EndpointRef,Int}
    next_connection_id::Int
end

Network() = Network(
    Dict{Symbol,AbstractComponent}(),
    Dict{Int,Connection}(),
    Dict{EndpointRef,Int}(),
    1,
)

function endpoint_exists(net::Network, ep::EndpointRef)
    comp = get(net.components, ep.component, nothing)
    comp === nothing && return false
    return has_port(comp, ep.port)
end

connection_id_for(net::Network, ep::EndpointRef) = get(net.connection_by_endpoint, ep, nothing)

function add_component!(net::Network, id::Symbol, component::AbstractComponent)
    haskey(net.components, id) && error("component id already exists: $id")
    net.components[id] = component
    return net
end

"""
Create a strict 2-port connection and return its connection id.
"""
function connect!(net::Network, a::EndpointRef, b::EndpointRef)
    endpoint_exists(net, a) || error("unknown endpoint: ($(a.component), $(a.port))")
    endpoint_exists(net, b) || error("unknown endpoint: ($(b.component), $(b.port))")
    a == b && error("connection endpoints must be distinct")

    connection_id_for(net, a) === nothing ||
        error("2-port rule violation: endpoint already connected ($(a.component), $(a.port))")
    connection_id_for(net, b) === nothing ||
        error("2-port rule violation: endpoint already connected ($(b.component), $(b.port))")

    cid = net.next_connection_id
    net.next_connection_id += 1

    conn = Connection(cid, a, b)
    net.connections_by_id[cid] = conn
    net.connection_by_endpoint[a] = cid
    net.connection_by_endpoint[b] = cid

    return cid
end

export EndpointRef, Connection
export Network, endpoint_exists, connection_id_for
export add_component!, connect!

end # module Network
