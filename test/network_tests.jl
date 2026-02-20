@testset "Network Structure" begin
    F = TurboMachineModel.Framework
    C = TurboMachineModel.Components
    P = TurboMachineModel.Physics.Fluids

    cmb = C.Combustor(0.04, 43e6, NamedTuple())
    shaft = C.InertialShaft(J=0.35, damping=0.01, n_ports=3)
    pm = P.PerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92],
    )
    tm = C.TurboMachineSection(mode=:compressor, performance_map=pm, eta_guess=0.9)

    net = F.Network()
    F.add_component!(net, :cmp, tm)
    F.add_component!(net, :cmb, cmb)
    F.add_component!(net, :shaft, shaft)
    F.connect!(net, F.EndpointRef(:cmp, :outlet), F.EndpointRef(:cmb, :inlet))
    F.connect!(net, F.EndpointRef(:cmp, :shaft), F.EndpointRef(:shaft, :shaft1))
    F.add_boundary!(
        net,
        :b_in,
        F.EndpointRef(:cmp, :inlet),
        :prescribed_state,
        (; pt=101_325.0, ht=300_000.0, mdot=-15.0, composition=:air),
    )
    F.add_boundary!(
        net,
        :b_out,
        F.EndpointRef(:cmb, :outlet),
        :prescribed_state,
        (; pt=99_000.0, ht=450_000.0, mdot=15.0, composition=:air),
    )
    F.add_boundary!(net, :b_shaft, F.EndpointRef(:shaft, :shaft2), :prescribed_speed, 10_000.0)
    F.add_boundary!(net, :b_shaft3, F.EndpointRef(:shaft, :shaft3), :prescribed_torque, 0.0)

    report = F.validate_network(net)
    @test report.valid

    adj = F.adjacency(net)
    @test length(adj[F.EndpointRef(:cmp, :outlet)]) == 1
    @test first(adj[F.EndpointRef(:cmp, :outlet)]) == F.EndpointRef(:cmb, :inlet)

    cmp_nbrs = F.component_neighbors(net, :cmp)
    @test :cmb in cmp_nbrs
    @test :shaft in cmp_nbrs

    fluid_conns = F.domain_connections(net, :fluid)
    shaft_conns = F.domain_connections(net, :shaft)
    @test length(fluid_conns) == 1
    @test length(shaft_conns) == 1

    @test_throws ErrorException F.connect!(
        net,
        F.EndpointRef(:cmp, :shaft),
        F.EndpointRef(:shaft, :shaft2),
    )
end
