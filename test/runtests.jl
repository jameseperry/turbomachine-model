using Test
using TurboMachineModel

@testset "TurboMachineModel.jl" begin
    @test isdefined(TurboMachineModel, :Framework)
    @test isdefined(TurboMachineModel, :Physics)
    @test isdefined(TurboMachineModel, :Components)
    @test isdefined(TurboMachineModel.Framework, :Model)
    @test isdefined(TurboMachineModel.Framework, :PortSpec)
    @test isdefined(TurboMachineModel.Components, :AbstractComponent)
    @test isdefined(TurboMachineModel.Components, :Combustor)
    @test isdefined(TurboMachineModel.Components, :TurboMachineSection)
    @test isdefined(TurboMachineModel.Components, :InertialShaft)

    cmb = TurboMachineModel.Components.Combustor(0.04, 43e6, NamedTuple())
    @test :inlet in keys(TurboMachineModel.Framework.port_specs(cmb))
    @test :outlet in keys(TurboMachineModel.Framework.port_specs(cmb))

    shaft = TurboMachineModel.Components.InertialShaft(
        J=0.35,
        damping=0.01,
        n_ports=3,
    )
    shaft_ports = keys(TurboMachineModel.Framework.port_specs(shaft))
    @test length(collect(shaft_ports)) == 3
    @test :shaft1 in shaft_ports
    @test :shaft2 in shaft_ports
    @test :shaft3 in shaft_ports

    pm = TurboMachineModel.Physics.Fluids.PerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],   # PR table
        [0.8 0.82; 0.9 0.92], # eta table
    )
    tm = TurboMachineModel.Components.TurboMachineSection(
        mode=:compressor,
        performance_map=pm,
        eta_guess=0.9,
    )
    @test tm.mode == :compressor
    @test :inlet in keys(TurboMachineModel.Framework.port_specs(tm))
    inlet_vars = first(TurboMachineModel.Framework.port_specs(tm)[:inlet].vars)
    @test inlet_vars.var == :Pt

    # Thermo utility smoke checks.
    Tt = 300.0
    Pt = 200_000.0
    M = 0.5
    gamma = 1.4
    T = TurboMachineModel.Physics.Fluids.static_temperature_from_total(Tt, M, gamma)
    P = TurboMachineModel.Physics.Fluids.static_pressure_from_total(Pt, M, gamma)
    @test isapprox(
        TurboMachineModel.Physics.Fluids.total_temperature_from_static(T, M, gamma),
        Tt;
        rtol=1e-12,
    )
    @test isapprox(
        TurboMachineModel.Physics.Fluids.total_pressure_from_static(P, M, gamma),
        Pt;
        rtol=1e-12,
    )

    # Performance map + correction smoke checks.
    vals = TurboMachineModel.Physics.Fluids.map_pr_eta(pm, 1.5, 15.0)
    @test isapprox(vals.PR, 3.5; rtol=1e-12)
    @test isapprox(vals.eta, 0.86; rtol=1e-12)

    from_inlet = TurboMachineModel.Physics.Fluids.map_pr_eta_from_stagnation(
        pm,
        10_000.0,
        15.0,
        300.0,
        100_000.0,
    )
    @test isapprox(from_inlet.N_corr, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.W_corr, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.PR, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    # Demo map helpers.
    cmp_demo = TurboMachineModel.Physics.Fluids.demo_compressor_map()
    trb_demo = TurboMachineModel.Physics.Fluids.demo_turbine_map()
    cmp_vals = TurboMachineModel.Physics.Fluids.map_pr_eta(cmp_demo, 0.8, 16.0)
    trb_vals = TurboMachineModel.Physics.Fluids.map_pr_eta(trb_demo, 0.8, 14.0)
    @test cmp_vals.PR > 1.0
    @test trb_vals.PR < 1.0
end
