@testset "Physics" begin
    P = TurboMachineModel.Physics.Fluids
    TP = TurboMachineModel.Physics
    TM = TP.Turbomachine.Compressor
    TT = TP.Turbomachine.Turbine
    AM = TP.Turbomachine.AxialMachine
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

    vals = TM.performance_from_stagnation(pm, 1.5, 15.0, 300.0, 100_000.0)
    @test isapprox(vals.PR, 3.5; rtol=1e-12)
    @test isapprox(vals.eta, 0.86; rtol=1e-12)
    @test vals.valid
    @test !vals.stall
    @test !vals.choke

    pm_limited = TM.TabulatedCompressorPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92];
        interpolation=:bilinear,
        mdot_corr_surge=[12.0, 14.0],
        mdot_corr_choke=[18.0, 19.0],
    )
    vals_stall = TM.performance_from_stagnation(pm_limited, 1.5, 12.0, 300.0, 100_000.0)
    vals_choke = TM.performance_from_stagnation(pm_limited, 1.5, 19.0, 300.0, 100_000.0)
    vals_mid = TM.performance_from_stagnation(pm_limited, 1.5, 15.5, 300.0, 100_000.0)
    @test vals_stall.stall
    @test !vals_stall.valid
    @test vals_choke.choke
    @test !vals_choke.valid
    @test vals_mid.valid
    pm_limited_domain = TM.performance_map_domain(pm_limited, 300.0, 100_000.0)
    @test isapprox(pm_limited_domain.mdot_flow_range.surge(1.5), 13.0; rtol=1e-12)
    @test isapprox(pm_limited_domain.mdot_flow_range.choke(1.5), 18.5; rtol=1e-12)

    pm_bicubic = TM.TabulatedCompressorPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
        [4.0 7.0 10.0; 6.0 9.0 12.0; 8.0 11.0 14.0],
        [0.70 0.80 0.90; 0.72 0.82 0.92; 0.74 0.84 0.94];
        interpolation=:bicubic,
    )
    vals_bicubic = TM.performance_from_stagnation(pm_bicubic, 2.5, 15.0, 300.0, 100_000.0)
    @test isapprox(vals_bicubic.PR, 8.5; rtol=1e-12)

    from_inlet = TM.performance_from_stagnation(pm, 10_000.0, 15.0, 300.0, 100_000.0)
    @test isapprox(from_inlet.speed_coord, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.flow_coord, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.PR, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    cmp_demo = TM.demo_tabulated_compressor_performance_map()
    cmp_vals = TM.performance_from_stagnation(
        cmp_demo,
        800.0,
        16.0,
        cmp_demo.Tt_ref,
        cmp_demo.Pt_ref,
    )
    @test cmp_vals.PR > 1.0
    cmp_domain = TM.performance_map_domain(cmp_demo, cmp_demo.Tt_ref, cmp_demo.Pt_ref)
    @test cmp_domain.omega == (600.0, 1000.0)
    @test cmp_domain.mdot == (12.0, 20.0)
    @test isapprox(cmp_domain.mdot_flow_range.surge(750.0), 12.0; rtol=1e-12)
    @test isapprox(cmp_domain.mdot_flow_range.choke(750.0), 20.0; rtol=1e-12)

    @testset "Compressor Map Coordinate Conversion" begin
        src = TM.demo_tabulated_compressor_performance_map(; interpolation=:bilinear)

        gamma = 1.4
        gas_constant = 287.05
        tip_radius_inlet = 0.22
        mean_radius_inlet = 0.18
        inlet_area = 0.060
        Tt_ref = src.Tt_ref
        Pt_ref = src.Pt_ref
        omega_ref_line = 800.0

        src_omega_grid = collect(U.table_xgrid(src.pr_map))
        src_mdot_corr_grid = collect(U.table_ygrid(src.pr_map))
        m_tip_grid = [
            omega * tip_radius_inlet / sqrt(gamma * gas_constant * Tt_ref) for omega in src_omega_grid
        ]
        rho0_ref = Pt_ref / (gas_constant * Tt_ref)
        phi_in_grid = [
            mdot / (rho0_ref * inlet_area * omega_ref_line * mean_radius_inlet) for mdot in src_mdot_corr_grid
        ]

        nd = TM.to_nondimensional_tabulated_compressor_map(
            src;
            gamma=gamma,
            gas_constant=gas_constant,
            tip_radius_inlet=tip_radius_inlet,
            mean_radius_inlet=mean_radius_inlet,
            inlet_area=inlet_area,
            Tt_in_ref=Tt_ref,
            Pt_in_ref=Pt_ref,
            m_tip_grid=m_tip_grid,
            phi_in_grid=phi_in_grid,
            interpolation=:bilinear,
        )
        @test nd isa TM.NonDimensionalTabulatedCompressorPerformanceMap

        src_vals = TM.performance_from_stagnation(src, 800.0, 16.0, Tt_ref, Pt_ref)
        nd_vals = TM.performance_from_stagnation(nd, 800.0, 16.0, Tt_ref, Pt_ref)
        @test isapprox(nd_vals.PR, src_vals.PR; rtol=1e-10)
        @test isapprox(nd_vals.eta, src_vals.eta; rtol=1e-10)

        back = TM.to_tabulated_compressor_map(
            nd;
            Tt_in_ref=Tt_ref,
            Pt_in_ref=Pt_ref,
            Tt_ref=Tt_ref,
            Pt_ref=Pt_ref,
            omega_corr_grid=src_omega_grid,
            mdot_corr_grid=src_mdot_corr_grid,
            interpolation=:bilinear,
        )
        @test back isa TM.TabulatedCompressorPerformanceMap

        back_vals = TM.performance_from_stagnation(back, 800.0, 16.0, Tt_ref, Pt_ref)
        src_grid_vals = TM.performance_from_stagnation(src, 800.0, 16.0, Tt_ref, Pt_ref)
        @test isapprox(back_vals.PR, src_grid_vals.PR; rtol=1e-10)
        @test isapprox(back_vals.eta, src_grid_vals.eta; rtol=1e-10)
        back_domain = TM.performance_map_domain(back, Tt_ref, Pt_ref)
        @test isapprox(back_domain.mdot_flow_range.surge(800.0), 12.0; rtol=1e-10)
        @test isapprox(back_domain.mdot_flow_range.choke(800.0), 20.0; rtol=1e-10)
    end

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
        domain = TT.performance_map_domain(map)
        @test domain.omega_corr == (1.0, 2.0)
        @test domain.pr_turb == (2.0, 3.0)

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
        table_map_toml_path = tempname() * ".toml"
        U.write_toml(pm.pr_map, table_map_toml_path)
        table_map_toml_loaded = U.read_toml(U.AbstractTableMap, table_map_toml_path)
        @test U.table_interpolation(table_map_toml_loaded) == :bilinear
        @test U.table_xgrid(table_map_toml_loaded) == U.table_xgrid(pm.pr_map)
        @test U.table_ygrid(table_map_toml_loaded) == U.table_ygrid(pm.pr_map)
        @test U.table_values(table_map_toml_loaded) == U.table_values(pm.pr_map)

        vals_ref = TM.performance_from_stagnation(pm_bicubic, 2.5, 15.0, 300.0, 100_000.0)
        compressor_map_toml_path = tempname() * ".toml"
        TM.write_toml(pm_bicubic, compressor_map_toml_path)
        pm_toml_loaded = TM.read_toml(TM.TabulatedCompressorPerformanceMap, compressor_map_toml_path)
        pm_generic_loaded = TM.read_performance_map_toml(compressor_map_toml_path; group="compressor_map")
        @test pm_generic_loaded isa TM.TabulatedCompressorPerformanceMap
        vals_toml_loaded = TM.performance_from_stagnation(pm_toml_loaded, 2.5, 15.0, 300.0, 100_000.0)
        @test isapprox(vals_toml_loaded.PR, vals_ref.PR; rtol=1e-12)
        @test isapprox(vals_toml_loaded.eta, vals_ref.eta; rtol=1e-12)

        pm_nd = TM.demo_nondimensional_tabulated_compressor_performance_map(; interpolation=:bicubic)
        compressor_map_nd_toml_path = tempname() * ".toml"
        TM.write_toml(pm_nd, compressor_map_nd_toml_path)
        pm_nd_toml_loaded = TM.read_toml(TM.NonDimensionalTabulatedCompressorPerformanceMap, compressor_map_nd_toml_path)
        pm_nd_generic_loaded = TM.read_performance_map_toml(compressor_map_nd_toml_path; group="compressor_map")
        @test pm_nd_generic_loaded isa TM.NonDimensionalTabulatedCompressorPerformanceMap
        vals_nd_ref = TM.performance_from_stagnation(pm_nd, 900.0, 12.0, 300.0, 101_325.0)
        vals_nd_loaded = TM.performance_from_stagnation(pm_nd_toml_loaded, 900.0, 12.0, 300.0, 101_325.0)
        @test isapprox(vals_nd_loaded.PR, vals_nd_ref.PR; rtol=1e-12)
        @test isapprox(vals_nd_loaded.eta, vals_nd_ref.eta; rtol=1e-12)

        meanline = TM.demo_compressor_meanline_model()
        m_mid = 0.5 * (meanline.m_tip_bounds[1] + meanline.m_tip_bounds[2])
        phi_mid = 0.5 * (meanline.phi_in_bounds[1] + meanline.phi_in_bounds[2])
        radii = AM.meanline_radii(meanline)
        @test length(radii) == length(meanline.rows)
        for (r, row) in zip(radii, meanline.rows)
            @test isapprox(r, 0.5 * (row.r_hub + row.r_tip); rtol=0, atol=1e-12)
        end
        @test AM.station_area(meanline, 1) ≈ π * (meanline.rows[1].r_tip^2 - meanline.rows[1].r_hub^2)
        @test AM.station_area(meanline, length(meanline.rows) + 1) ≈ π * (meanline.rows[end].r_tip^2 - meanline.rows[end].r_hub^2)
        if length(meanline.rows) >= 2
            expected_a2 = 0.5 * (
                π * (meanline.rows[1].r_tip^2 - meanline.rows[1].r_hub^2) +
                π * (meanline.rows[2].r_tip^2 - meanline.rows[2].r_hub^2)
            )
            @test AM.station_area(meanline, 2) ≈ expected_a2
        end
        meanline_vals = AM.streamtube_solve_with_phi(meanline, m_mid, phi_mid)
        @test isfinite(meanline_vals.PR)
        @test isfinite(meanline_vals.eta)
        first_rotor = findfirst(row -> row.kind == :rotor, meanline.rows)
        idx_ref = isnothing(first_rotor) ? 1 : first_rotor
        nu_u_ref = meanline.rows[idx_ref].speed_ratio_to_ref * m_mid * radii[idx_ref] / meanline.r_tip_ref
        nu_x_mid = phi_mid * abs(nu_u_ref)
        meanline_vals_core = AM.streamtube_solve(meanline, radii, m_mid, nu_x_mid, 0.0)
        @test isfinite(meanline_vals_core.PR)
        @test isfinite(meanline_vals_core.eta)
        @test isapprox(meanline_vals_core.PR, meanline_vals.PR; rtol=1e-12)
        @test isapprox(meanline_vals_core.eta, meanline_vals.eta; rtol=1e-12)

        nd_from_meanline = TM.tabulate_compressor_meanline_model(
            meanline;
            n_speed=9,
            n_flow=11,
            interpolation=:bilinear,
        )
        @test nd_from_meanline isa TM.NonDimensionalTabulatedCompressorPerformanceMap
        m_grid = U.table_xgrid(nd_from_meanline.pr_map)
        phi_grid = U.table_ygrid(nd_from_meanline.pr_map)
        m_sample = m_grid[4]
        phi_sample = phi_grid[6]
        r_tip_1 = meanline.r_tip_ref
        r_mean_1 = radii[idx_ref]
        a0 = sqrt(meanline.gamma * meanline.gas_constant * 300.0)
        rho0 = 101_325.0 / (meanline.gas_constant * 300.0)
        A_phys_1 = AM.station_area(meanline, 1)
        omega_sample = m_sample * a0 / r_tip_1
        mdot_sample = phi_sample * abs(omega_sample) * r_mean_1 * rho0 * A_phys_1
        vals_direct = AM.streamtube_solve_with_phi(
            meanline,
            m_sample,
            phi_sample,
        )
        vals_nd_table = TM.performance_from_stagnation(
            nd_from_meanline,
            omega_sample,
            mdot_sample,
            300.0,
            101_325.0,
        )
        @test isapprox(vals_nd_table.PR, vals_direct.PR; rtol=1e-9)
        @test isapprox(vals_nd_table.eta, vals_direct.eta; rtol=1e-9)

        meanline_path = tempname() * ".toml"
        TM.write_toml(meanline, meanline_path)
        meanline_loaded = TM.read_toml(TM.CompressorMeanlineModel, meanline_path)
        @test_throws ErrorException TM.read_performance_map_toml(meanline_path; group="compressor_meanline_model")
        meanline_vals_loaded = AM.streamtube_solve_with_phi(meanline_loaded, m_mid, phi_mid)
        @test isapprox(meanline_vals_loaded.PR, meanline_vals.PR; rtol=1e-12)
        @test isapprox(meanline_vals_loaded.eta, meanline_vals.eta; rtol=1e-12)

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
