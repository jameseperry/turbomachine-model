@testset "Network Structure" begin
    F = TurboMachineModel.Network
    C = TurboMachineModel.Components
    P = TurboMachineModel.Physics.Fluids

    cmb = C.Combustor(0.04, 43e6, NamedTuple())
    shaft = C.InertialShaft(J=0.35, damping=0.01, n_ports=1)
    pm = P.PerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92],
    )
    tm = C.Turbomachine(mode=:compressor, performance_map=pm, eta_guess=0.9)

    net = F.Network()
    F.add_component!(net, :cmp, tm)
    F.add_component!(net, :cmb, cmb)
    F.add_component!(net, :shaft, shaft)

    c_fluid_1 = F.connect!(net, F.EndpointRef(:cmp, :outlet), F.EndpointRef(:cmb, :inlet))
    c_fluid_2 = F.connect!(net, F.EndpointRef(:cmb, :outlet), F.EndpointRef(:cmp, :inlet))
    c_shaft = F.connect!(net, F.EndpointRef(:cmp, :shaft), F.EndpointRef(:shaft, :shaft1))

    @test c_fluid_1 == 1
    @test c_fluid_2 == 2
    @test c_shaft == 3
    @test length(net.connections_by_id) == 3

    @test_throws ErrorException F.connect!(
        net,
        F.EndpointRef(:cmp, :shaft),
        F.EndpointRef(:shaft, :shaft1),
    )
end
