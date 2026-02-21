@testset "Physics" begin
    P = TurboMachineModel.Physics.Fluids
    TP = TurboMachineModel.Physics

    Tt = 300.0
    Pt = 200_000.0
    M = 0.5
    gamma = 1.4
    T = P.static_temperature_from_total(Tt, M, gamma)
    Pstat = P.static_pressure_from_total(Pt, M, gamma)
    @test isapprox(P.total_temperature_from_static(T, M, gamma), Tt; rtol=1e-12)
    @test isapprox(P.total_pressure_from_static(Pstat, M, gamma), Pt; rtol=1e-12)

    # Velocity from (p, h, mdot, A) through density closure.
    rho_fn = (p, h) -> 2.0 + 1e-8 * p + 0.0 * h
    v1 = P.velocity_from_ph_mdot(100_000.0, 300_000.0, 10.0, 0.5, rho_fn)
    @test isapprox(v1, 10.0 / ((2.0 + 1e-8 * 100_000.0) * 0.5); rtol=1e-12)
    @test isapprox(P.velocity_from_massflow(10.0, 2.001, 0.5), v1; rtol=1e-12)
    @test_throws ErrorException P.velocity_from_massflow(1.0, 0.0, 0.5)
    @test_throws ErrorException P.velocity_from_massflow(1.0, 1.0, 0.0)

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
    @test isapprox(from_inlet.omega_corr, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.mdot_corr, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.PR, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    cmp_demo = P.demo_compressor_map()
    trb_demo = P.demo_turbine_map()
    cmp_vals = P.map_pr_eta(cmp_demo, 0.8, 16.0)
    trb_vals = P.map_pr_eta(trb_demo, 0.8, 14.0)
    @test cmp_vals.PR > 1.0
    @test trb_vals.PR < 1.0

    @testset "Turbomachine Residuals" begin
        eos = P.ideal_EOS()[:air]
        map = P.PerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [10.0, 20.0],
            [2.0 2.0; 2.0 2.0],
            [0.8 0.8; 0.8 0.8],
        )
        pt_in = 100_000.0
        ht_in = 300_000.0
        mdot = 15.0
        omega = 12_000.0
        pt_out = 2.0 * pt_in
        h2s = P.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
        ht_out = ht_in + (h2s - ht_in) / 0.8
        tau = mdot * (ht_out - ht_in) / omega

        R_pr, R_eta, R_P = TP.turbomachine_residuals(
            map,
            eos,
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau,
        )
        @test isapprox(R_pr, 0.0; atol=1e-8)
        @test isapprox(R_eta, 0.0; atol=1e-8)
        @test isapprox(R_P, 0.0; atol=1e-8)
    end
end
