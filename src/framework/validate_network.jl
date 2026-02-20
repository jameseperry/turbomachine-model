function _issue!(issues::Vector{ValidationIssue}, level::Symbol, code::Symbol, msg::String)
    push!(issues, ValidationIssue(level, code, msg))
end

_is_error(i::ValidationIssue) = i.level == :error

function _port_spec(net::Network, ep::EndpointRef)
    comp = get(net.components, ep.component, nothing)
    comp === nothing && return nothing
    return get(port_specs(comp), ep.port, nothing)
end

function _varspec_dict(vars)
    Dict(v.var => (unit=v.unit, coupling=v.coupling) for v in vars)
end

function _validate_connection_vars!(
    issues::Vector{ValidationIssue},
    conn::ConnectionSpec,
    from_spec,
    to_spec,
)
    from_vars = _varspec_dict(from_spec.vars)
    to_vars = _varspec_dict(to_spec.vars)

    used_vars = conn.vars
    if conn.vars == :all
        if Set(keys(from_vars)) != Set(keys(to_vars))
            _issue!(
                issues,
                :error,
                :vars_mismatch_all,
                "connection :all requires matching variable sets ($(conn.from.component).$(conn.from.port) <-> $(conn.to.component).$(conn.to.port))",
            )
            return
        end
        used_vars = collect(keys(from_vars))
    end

    for v in used_vars
        if !haskey(from_vars, v)
            _issue!(
                issues,
                :error,
                :missing_var_source,
                "missing variable :$v on source endpoint ($(conn.from.component), $(conn.from.port))",
            )
            continue
        end
        if !haskey(to_vars, v)
            _issue!(
                issues,
                :error,
                :missing_var_target,
                "missing variable :$v on target endpoint ($(conn.to.component), $(conn.to.port))",
            )
            continue
        end

        if from_vars[v].coupling != to_vars[v].coupling
            _issue!(issues, :error, :coupling_mismatch, "coupling mismatch for variable :$v")
        end
        if from_vars[v].coupling ∉ (:equality, :conservation)
            _issue!(issues, :error, :invalid_coupling, "invalid coupling for variable :$v")
        end
        if from_vars[v].unit != to_vars[v].unit
            _issue!(issues, :error, :unit_mismatch, "unit mismatch for variable :$v")
        end
    end
end

function _validate_connection!(
    issues::Vector{ValidationIssue},
    net::Network,
    conn::ConnectionSpec,
)
    from_spec = _port_spec(net, conn.from)
    to_spec = _port_spec(net, conn.to)

    from_spec === nothing &&
        _issue!(issues, :error, :unknown_endpoint, "unknown connection endpoint: ($(conn.from.component), $(conn.from.port))")
    to_spec === nothing &&
        _issue!(issues, :error, :unknown_endpoint, "unknown connection endpoint: ($(conn.to.component), $(conn.to.port))")

    (from_spec === nothing || to_spec === nothing) && return

    from_spec.domain == to_spec.domain ||
        _issue!(issues, :error, :domain_mismatch, "domain mismatch for connection ($(conn.from.component).$(conn.from.port) <-> $(conn.to.component).$(conn.to.port))")

    conn.from != conn.to ||
        _issue!(issues, :error, :self_connection, "self-connection is not allowed")

    _validate_connection_vars!(issues, conn, from_spec, to_spec)
end

function _validate_2port_rule!(issues::Vector{ValidationIssue}, net::Network)
    endpoint_counts = Dict{EndpointRef,Int}()
    for c in net.connections
        endpoint_counts[c.from] = get(endpoint_counts, c.from, 0) + 1
        endpoint_counts[c.to] = get(endpoint_counts, c.to, 0) + 1
    end
    for (ep, n) in endpoint_counts
        if n > 1
            _issue!(
                issues,
                :error,
                :two_port_rule,
                "2-port rule violation at endpoint ($(ep.component), $(ep.port)); found $n connections",
            )
        end
    end
end

function _validate_required_ports!(issues::Vector{ValidationIssue}, net::Network)
    connected = Set{EndpointRef}()
    for c in net.connections
        push!(connected, c.from)
        push!(connected, c.to)
    end
    bounded = Set{EndpointRef}(b.target for b in net.boundaries)

    for (id, comp) in net.components
        specs = port_specs(comp)
        for p in required_ports(comp)
            if !haskey(specs, p)
                _issue!(issues, :error, :missing_required_port_spec, "component $id requires unknown port :$p")
                continue
            end
            ep = EndpointRef(id, p)
            if ep ∉ connected && ep ∉ bounded
                _issue!(issues, :error, :unconstrained_required_port, "required port not connected/bounded: ($id, $p)")
            end
        end
    end
end

function _validate_boundaries!(issues::Vector{ValidationIssue}, net::Network)
    ids = Set{Symbol}()
    target_types = Set{Tuple{EndpointRef,Symbol}}()
    for b in net.boundaries
        if b.id in ids
            _issue!(issues, :error, :duplicate_boundary_id, "duplicate boundary id: $(b.id)")
        else
            push!(ids, b.id)
        end

        _port_spec(net, b.target) === nothing &&
            _issue!(
                issues,
                :error,
                :unknown_boundary_target,
                "unknown boundary target endpoint: ($(b.target.component), $(b.target.port))",
            )

        key = (b.target, b.bc_type)
        if key in target_types
            _issue!(
                issues,
                :error,
                :duplicate_boundary_target,
                "duplicate boundary for same target/bc_type: ($(b.target.component), $(b.target.port)) / $(b.bc_type)",
            )
        else
            push!(target_types, key)
        end
    end
end

function _validate_component_configs!(issues::Vector{ValidationIssue}, net::Network)
    for (id, comp) in net.components
        try
            validate(comp)
        catch e
            _issue!(issues, :error, :component_invalid, "component $id failed validate(): $(sprint(showerror, e))")
        end
    end
end

function _warn_disconnected_islands!(issues::Vector{ValidationIssue}, net::Network)
    ids = collect(keys(net.components))
    isempty(ids) && return

    nbrs = Dict(id => Set{Symbol}() for id in ids)
    for c in net.connections
        push!(nbrs[c.from.component], c.to.component)
        push!(nbrs[c.to.component], c.from.component)
    end

    seen = Set{Symbol}()
    n_islands = 0
    for id in ids
        id in seen && continue
        n_islands += 1
        stack = [id]
        while !isempty(stack)
            cur = pop!(stack)
            cur in seen && continue
            push!(seen, cur)
            append!(stack, collect(nbrs[cur]))
        end
    end

    if n_islands > 1
        _issue!(
            issues,
            :warning,
            :disconnected_islands,
            "network has $n_islands disconnected component islands",
        )
    end
end

function validate_network(net::Network)
    issues = ValidationIssue[]

    _validate_component_configs!(issues, net)
    for c in net.connections
        _validate_connection!(issues, net, c)
    end
    _validate_2port_rule!(issues, net)
    _validate_boundaries!(issues, net)
    _validate_required_ports!(issues, net)
    _warn_disconnected_islands!(issues, net)

    return ValidationReport(!any(_is_error, issues), issues)
end
