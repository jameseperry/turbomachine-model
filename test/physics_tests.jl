@testset "Physics" begin
    P = TurboMachineModel.Physics.Fluids
    TP = TurboMachineModel.Physics
    TM = TP.Turbomachine

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

    pm = TM.TabulatedPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92],
    )

    vals = TM.performance_map(pm, 1.5, 15.0)
    @test isapprox(vals.PR, 3.5; rtol=1e-12)
    @test isapprox(vals.eta, 0.86; rtol=1e-12)

    from_inlet = TM.performance_map_from_stagnation(pm, 10_000.0, 15.0, 300.0, 100_000.0)
    @test isapprox(from_inlet.omega_corr, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.mdot_corr, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.PR, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    cmp_demo = TM.demo_compressor_performance_map()
    trb_demo = TM.demo_turbine_performance_map()
    cmp_vals = TM.performance_map(cmp_demo, 0.8, 16.0)
    trb_vals = TM.performance_map(trb_demo, 0.8, 14.0)
    @test cmp_vals.PR > 1.0
    @test trb_vals.PR < 1.0

    @testset "Turbomachine Residuals" begin
        eos = P.ideal_EOS()[:air]
        map = TM.TabulatedPerformanceMap(
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

        R_pout, R_dh_eff, R_P = TM.turbomachine_residuals(
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
        @test isapprox(R_pout, 0.0; atol=1e-8)
        @test isapprox(R_dh_eff, 0.0; atol=1e-8)
        @test isapprox(R_P, 0.0; atol=1e-8)
    end

    @testset "Turbomachine Scaled Residuals" begin
        eos = P.ideal_EOS()[:air]
        map = TM.TabulatedPerformanceMap(
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
        pt_out = 210_000.0
        ht_out = 380_000.0
        tau = 120.0

        R_pout, R_dh_eff, R_P = TM.turbomachine_residuals(
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

        scales = TM.turbomachine_residual_scales(
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau,
        )
        r_pout, r_dh_eff, r_P = TM.turbomachine_residuals_scaled(
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
        @test isapprox(r_pout, R_pout / scales.pressure_scale; rtol=1e-12)
        @test isapprox(r_dh_eff, R_dh_eff / scales.enthalpy_scale; rtol=1e-12)
        @test isapprox(r_P, R_P / scales.power_scale; rtol=1e-12)

        r2_pout, r2_dh_eff, r2_P = TM.turbomachine_residuals_scaled(
            map,
            eos,
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau;
            pressure_scale=2.0e5,
            enthalpy_scale=5.0e5,
            power_scale=3.0e6,
        )
        @test isapprox(r2_pout, R_pout / 2.0e5; rtol=1e-12)
        @test isapprox(r2_dh_eff, R_dh_eff / 5.0e5; rtol=1e-12)
        @test isapprox(r2_P, R_P / 3.0e6; rtol=1e-12)
    end

    @testset "Turbomachine Operating Point Solve" begin
        eos = P.ideal_EOS()[:air]
        map = TM.TabulatedPerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [10.0, 20.0],
            [1.8 2.2; 1.8 2.2],
            [0.8 0.8; 0.8 0.8],
        )

        pt_in = 100_000.0
        ht_in = 300_000.0
        pt_out = 200_000.0
        omega = 12_000.0

        Tt_in = P.temperature(eos, pt_in, ht_in)
        mdot_expected = 15.0 * (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
        h2s = P.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
        ht_out_expected = ht_in + (h2s - ht_in) / 0.8
        tau_expected = mdot_expected * (ht_out_expected - ht_in) / omega

        sol = TM.solve_turbomachine_operating_point(
            map,
            eos;
            pt_in=pt_in,
            ht_in=ht_in,
            pt_out=pt_out,
            omega=omega,
            mdot_guess=mdot_expected * 0.9,
            ht_out_guess=ht_out_expected * 1.05,
            tau_guess=tau_expected * 0.8,
        )

        @test sol.converged
        @test isapprox(sol.mdot, mdot_expected; rtol=1e-8)
        @test isapprox(sol.ht_out, ht_out_expected; rtol=1e-8)
        @test isapprox(sol.tau, tau_expected; rtol=1e-8)
        @test isapprox(sol.residuals[1], 0.0; atol=1e-9)
        @test isapprox(sol.residuals[2], 0.0; atol=1e-9)
        @test isapprox(sol.residuals[3], 0.0; atol=1e-9)
    end
end
