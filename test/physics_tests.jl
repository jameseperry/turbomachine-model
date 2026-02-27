@testset "Physics" begin
    P = TurboMachineModel.Physics.Fluids
    TP = TurboMachineModel.Physics
    TM = TP.Turbomachine.Compressor
    TT = TP.Turbomachine.Turbine
    U = TurboMachineModel.Utility

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

    pm = TM.TabulatedCompressorPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92];
        interpolation=:bilinear,
    )

    vals = TM.compressor_performance_map(pm, 1.5, 15.0)
    @test isapprox(vals.PR, 3.5; rtol=1e-12)
    @test isapprox(vals.eta, 0.86; rtol=1e-12)

    pm_bicubic = TM.TabulatedCompressorPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
        [4.0 7.0 10.0; 6.0 9.0 12.0; 8.0 11.0 14.0],
        [0.70 0.80 0.90; 0.72 0.82 0.92; 0.74 0.84 0.94];
        interpolation=:bicubic,
    )
    vals_bicubic = TM.compressor_performance_map(pm_bicubic, 2.5, 15.0)
    @test isapprox(vals_bicubic.PR, 8.5; rtol=1e-12)

    from_inlet = TM.compressor_performance_map_from_stagnation(pm, 10_000.0, 15.0, 300.0, 100_000.0)
    @test isapprox(from_inlet.omega_corr, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.mdot_corr, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.PR, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    cmp_demo = TM.demo_compressor_performance_map()
    cmp_vals = TM.compressor_performance_map(cmp_demo, 0.8, 16.0)
    @test cmp_vals.PR > 1.0

    @testset "Turbomachine Residuals" begin
        eos = P.ideal_EOS()[:air]
        map = TM.TabulatedCompressorPerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [10.0, 20.0],
            [2.0 2.0; 2.0 2.0],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
        )
        pt_in = 100_000.0
        ht_in = 300_000.0
        mdot = 15.0
        omega = 12_000.0
        pt_out = 2.0 * pt_in
        h2s = P.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
        ht_out = ht_in + (h2s - ht_in) / 0.8
        tau = mdot * (ht_out - ht_in) / omega

        R_pout, R_dh_eff, R_P = TM.compressor_residuals(
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
        map = TM.TabulatedCompressorPerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [10.0, 20.0],
            [2.0 2.0; 2.0 2.0],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
        )
        pt_in = 100_000.0
        ht_in = 300_000.0
        mdot = 15.0
        omega = 12_000.0
        pt_out = 210_000.0
        ht_out = 380_000.0
        tau = 120.0

        R_pout, R_dh_eff, R_P = TM.compressor_residuals(
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

        scales = TM.compressor_residual_scales(
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau,
        )
        r_pout, r_dh_eff, r_P = TM.compressor_residuals_scaled(
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

        r2_pout, r2_dh_eff, r2_P = TM.compressor_residuals_scaled(
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
        map = TM.TabulatedCompressorPerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [10.0, 20.0],
            [1.8 2.2; 1.8 2.2],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
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

        sol = TM.solve_compressor_operating_point(
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

    @testset "Turbine APIs" begin
        eos = P.ideal_EOS()[:air]
        map = TT.TabulatedTurbinePerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [2.0, 3.0],
            [10.0 12.0; 14.0 16.0],
            [0.85 0.86; 0.87 0.88];
            interpolation=:bilinear,
        )

        vals = TT.turbine_performance_map(map, 1.5, 2.5)
        @test isapprox(vals.mdot_corr, 13.0; rtol=1e-12)
        @test isapprox(vals.eta, 0.865; rtol=1e-12)

        from_stag = TT.turbine_performance_map_from_stagnation(map, 1.5, 200_000.0, 100_000.0, 300.0)
        @test isapprox(from_stag.omega_corr, 1.5; rtol=1e-12)
        @test isapprox(from_stag.PR_turb, 2.0; rtol=1e-12)
        @test isapprox(from_stag.mdot_corr, 12.0; rtol=1e-12)
        @test isapprox(from_stag.mdot, 24.0; rtol=1e-12)
        @test isapprox(from_stag.eta, 0.86; rtol=1e-12)

        map_const = TT.TabulatedTurbinePerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [2.0, 3.0],
            [10.0 10.0; 10.0 10.0],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
        )

        pt_in = 100_000.0
        ht_in = 300_000.0
        pt_out = 50_000.0
        omega = 1.5
        Tt_in = P.temperature(eos, pt_in, ht_in)
        map_vals = TT.turbine_performance_map_from_stagnation(map_const, omega, pt_in, pt_out, Tt_in)
        mdot = map_vals.mdot
        h2s = P.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
        ht_out = ht_in - map_vals.eta * (ht_in - h2s)
        tau = mdot * (ht_out - ht_in) / omega

        R_mdot_map, R_dh_eff, R_P = TT.turbine_residuals(
            map_const,
            eos,
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau,
        )
        @test isapprox(R_mdot_map, 0.0; atol=1e-8)
        @test isapprox(R_dh_eff, 0.0; atol=1e-8)
        @test isapprox(R_P, 0.0; atol=1e-8)

        sol = TT.solve_turbine_operating_point(
            map_const,
            eos;
            pt_in=pt_in,
            ht_in=ht_in,
            pt_out=pt_out,
            omega=omega,
            mdot_guess=mdot * 1.1,
            ht_out_guess=ht_out * 0.95,
            tau_guess=tau * 1.2,
        )
        @test sol.converged
        @test isapprox(sol.mdot, mdot; rtol=1e-8)
        @test isapprox(sol.ht_out, ht_out; rtol=1e-8)
        @test isapprox(sol.tau, tau; rtol=1e-8)
    end

    @testset "Map IO" begin
        analytic_map = TM.AnalyticCompressorPerformanceMap{Float64}(
            Pi_max=1.75,
            pr_speed_exp=2.25,
            eta_max=0.90,
            Tt_ref=300.0,
            Pt_ref=95_000.0,
        )
        analytic_map_path = tempname() * ".toml"
        TM.write_toml(analytic_map, analytic_map_path)
        analytic_map_loaded = TM.read_toml(TM.AnalyticCompressorPerformanceMap, analytic_map_path)
        for field in fieldnames(TM.AnalyticCompressorPerformanceMap{Float64})
            @test isapprox(getfield(analytic_map_loaded, field), getfield(analytic_map, field); rtol=1e-12)
        end

        spec = TM.CompressorSpec(
            pr_design=4.0,
            eta_design=0.86,
            flow_range=0.45,
            surge_margin=0.62,
            choke_sharpness=0.35,
            speed_sensitivity=0.70,
        )
        spec_path = tempname() * ".toml"
        TM.write_toml(spec, spec_path)
        spec_loaded = TM.read_toml(TM.CompressorSpec, spec_path)
        for field in fieldnames(TM.CompressorSpec{Float64})
            @test isapprox(getfield(spec_loaded, field), getfield(spec, field); rtol=1e-12)
        end

        design = TM.CompressorDesign(
            kind=:centrifugal,
            stage_count=1,
            stage_loading=0.65,
            tip_mach_design=0.75,
            diffusion_aggressiveness=0.45,
            clearance_fraction=0.15,
            diffuser_quality=0.72,
            variable_geometry=0.10,
            reynolds_quality=0.88,
        )
        design_path = tempname() * ".toml"
        TM.write_toml(design, design_path)
        design_loaded = TM.read_toml(TM.CompressorDesign, design_path)
        @test design_loaded.kind == design.kind
        @test design_loaded.stage_count == design.stage_count
        for field in fieldnames(TM.CompressorDesign{Float64})
            field in (:kind, :stage_count) && continue
            @test isapprox(getfield(design_loaded, field), getfield(design, field); rtol=1e-12)
        end

        table_map_toml_path = tempname() * ".toml"
        U.write_toml(pm.pr_map, table_map_toml_path)
        table_map_toml_loaded = U.read_toml(U.AbstractTableMap, table_map_toml_path)
        @test U.table_interpolation(table_map_toml_loaded) == :bilinear
        @test U.table_xgrid(table_map_toml_loaded) == U.table_xgrid(pm.pr_map)
        @test U.table_ygrid(table_map_toml_loaded) == U.table_ygrid(pm.pr_map)
        @test U.table_values(table_map_toml_loaded) == U.table_values(pm.pr_map)

        vals_ref = TM.compressor_performance_map(pm_bicubic, 2.5, 15.0)
        compressor_map_toml_path = tempname() * ".toml"
        TM.write_toml(pm_bicubic, compressor_map_toml_path)
        pm_toml_loaded = TM.read_toml(TM.TabulatedCompressorPerformanceMap, compressor_map_toml_path)
        vals_toml_loaded = TM.compressor_performance_map(pm_toml_loaded, 2.5, 15.0)
        @test isapprox(vals_toml_loaded.PR, vals_ref.PR; rtol=1e-12)
        @test isapprox(vals_toml_loaded.eta, vals_ref.eta; rtol=1e-12)

        turbine_map = TT.TabulatedTurbinePerformanceMap(
            300.0,
            100_000.0,
            [1.0, 2.0],
            [2.0, 3.0],
            [10.0 10.0; 10.0 10.0],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
        )
        vals_t_ref = TT.turbine_performance_map(turbine_map, 1.5, 2.5)
        turbine_map_toml_path = tempname() * ".toml"
        TT.write_toml(turbine_map, turbine_map_toml_path)
        turbine_map_toml_loaded = TT.read_toml(TT.TabulatedTurbinePerformanceMap, turbine_map_toml_path)
        vals_t_toml_loaded = TT.turbine_performance_map(turbine_map_toml_loaded, 1.5, 2.5)
        @test isapprox(vals_t_toml_loaded.mdot_corr, vals_t_ref.mdot_corr; rtol=1e-12)
        @test isapprox(vals_t_toml_loaded.eta, vals_t_ref.eta; rtol=1e-12)
    end
end
