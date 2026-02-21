@testset "Components" begin
    F = TurboMachineModel.Framework
    C = TurboMachineModel.Components
    P = TurboMachineModel.Physics.Fluids

    cmb = C.Combustor(0.04, 43e6, NamedTuple())
    @test :inlet in keys(F.port_specs(cmb))
    @test :outlet in keys(F.port_specs(cmb))

    pln = C.Plenum(volume=0.25)
    @test :inlet in keys(F.port_specs(pln))
    @test :outlet in keys(F.port_specs(pln))
    @test F.port_specs(pln)[:inlet].domain == :fluid
    @test F.port_specs(pln)[:outlet].domain == :fluid
    @test_throws ErrorException C.Plenum(volume=0.0)

    shaft = C.InertialShaft(J=0.35, damping=0.01, n_ports=3)
    shaft_ports = keys(F.port_specs(shaft))
    @test length(collect(shaft_ports)) == 3
    @test :shaft1 in shaft_ports
    @test :shaft2 in shaft_ports
    @test :shaft3 in shaft_ports

    gb = C.Gearbox(ratio=2.5, efficiency=0.98)
    gb_ports = keys(F.port_specs(gb))
    @test :input in gb_ports
    @test :output in gb_ports
    @test F.port_specs(gb)[:input].domain == :shaft
    @test F.port_specs(gb)[:output].domain == :shaft
    @test_throws ErrorException C.Gearbox(ratio=0.0)
    @test_throws ErrorException C.Gearbox(ratio=2.0, efficiency=0.0)

    pm = P.PerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92],
    )
    tm = C.TurboMachineSection(
        mode=:compressor,
        performance_map=pm,
        eta_guess=0.9,
    )
    @test tm.mode == :compressor
    @test :inlet in keys(F.port_specs(tm))
    inlet_vars = first(F.port_specs(tm)[:inlet].vars)
    @test inlet_vars.var == :pt
end
