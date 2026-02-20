function adjacency(net::Network)
    adj = Dict{EndpointRef,Vector{EndpointRef}}()
    for c in net.connections
        push!(get!(adj, c.from, EndpointRef[]), c.to)
        push!(get!(adj, c.to, EndpointRef[]), c.from)
    end
    return adj
end

function connections_for(net::Network, endpoint::EndpointRef)
    filter(c -> c.from == endpoint || c.to == endpoint, net.connections)
end

function component_neighbors(net::Network, component_id::Symbol)
    nbrs = Set{Symbol}()
    for c in net.connections
        if c.from.component == component_id
            push!(nbrs, c.to.component)
        elseif c.to.component == component_id
            push!(nbrs, c.from.component)
        end
    end
    return collect(nbrs)
end

function _port_spec_for(net::Network, ep::EndpointRef)
    comp = net.components[ep.component]
    return port_specs(comp)[ep.port]
end

function domain_connections(net::Network, domain::Symbol)
    filter(net.connections) do c
        _port_spec_for(net, c.from).domain == domain &&
            _port_spec_for(net, c.to).domain == domain
    end
end
