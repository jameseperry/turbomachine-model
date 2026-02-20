@testset "Components" begin
    F = TurboMachineModel.Framework
    C = TurboMachineModel.Components
    P = TurboMachineModel.Physics.Fluids

    cmb = C.Combustor(0.04, 43e6, NamedTuple())
    @test :inlet in keys(F.port_specs(cmb))
    @test :outlet in keys(F.port_specs(cmb))

    shaft = C.InertialShaft(J=0.35, damping=0.01, n_ports=3)
    shaft_ports = keys(F.port_specs(shaft))
    @test length(collect(shaft_ports)) == 3
    @test :shaft1 in shaft_ports
    @test :shaft2 in shaft_ports
    @test :shaft3 in shaft_ports

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
    @test inlet_vars.var == :p
end
