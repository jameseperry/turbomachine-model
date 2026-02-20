@testset "Physics" begin
    P = TurboMachineModel.Physics.Fluids

    Tt = 300.0
    Pt = 200_000.0
    M = 0.5
    gamma = 1.4
    T = P.static_temperature_from_total(Tt, M, gamma)
    Pstat = P.static_pressure_from_total(Pt, M, gamma)
    @test isapprox(P.total_temperature_from_static(T, M, gamma), Tt; rtol=1e-12)
    @test isapprox(P.total_pressure_from_static(Pstat, M, gamma), Pt; rtol=1e-12)

    pm = P.PerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92],
    )

    vals = P.map_pr_eta(pm, 1.5, 15.0)
    @test isapprox(vals.PR, 3.5; rtol=1e-12)
    @test isapprox(vals.eta, 0.86; rtol=1e-12)

    from_inlet = P.map_pr_eta_from_stagnation(pm, 10_000.0, 15.0, 300.0, 100_000.0)
    @test isapprox(from_inlet.N_corr, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.W_corr, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.PR, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    cmp_demo = P.demo_compressor_map()
    trb_demo = P.demo_turbine_map()
    cmp_vals = P.map_pr_eta(cmp_demo, 0.8, 16.0)
    trb_vals = P.map_pr_eta(trb_demo, 0.8, 14.0)
    @test cmp_vals.PR > 1.0
    @test trb_vals.PR < 1.0
end
