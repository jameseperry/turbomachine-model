@testset "Physics" begin
    P = TurboMachineModel.Physics.Fluids
    TP = TurboMachineModel.Physics
    TM = TP.Turbomachine
    TTM = TP.Turbomachine
    AM = TP.Turbomachine.AxialMachine
    U = TurboMachineModel.Utility

    inlet_Vx(model, eos, pt_in, ht_in, mdot, Vtheta_inlet) = begin
        inlet = P.velocity_from_stagnation_massflow(
            eos,
            mdot,
            pt_in,
            ht_in,
            AM.station_area(model, 1),
            Vtheta_inlet;
            prefer=:low,
        )
        @test inlet.converged
        inlet.Vx
    end

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
        [0.8 0.82; 0.9 0.92];
        interpolation=:bilinear,
    )

    vals = TM.performance_from_stagnation(pm, 1.5, 15.0, 300.0, 100_000.0)
    @test isapprox(vals.pressure_ratio, 3.5; rtol=1e-12)
    @test isapprox(vals.eta, 0.86; rtol=1e-12)
    @test vals.valid
    @test !vals.low_flow
    @test !vals.high_flow

    pm_limited = TM.TabulatedPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92];
        interpolation=:bilinear,
        flow_min=[12.0, 14.0],
        flow_max=[18.0, 19.0],
    )
    vals_stall = TM.performance_from_stagnation(pm_limited, 1.5, 12.0, 300.0, 100_000.0)
    vals_choke = TM.performance_from_stagnation(pm_limited, 1.5, 19.0, 300.0, 100_000.0)
    vals_mid = TM.performance_from_stagnation(pm_limited, 1.5, 15.5, 300.0, 100_000.0)
    @test vals_stall.low_flow
    @test !vals_stall.valid
    @test vals_choke.high_flow
    @test !vals_choke.valid
    @test vals_mid.valid
    pm_limited_domain = TM.performance_map_domain(pm_limited, 300.0, 100_000.0)
    @test isapprox(pm_limited_domain.mdot_flow_range.min(1.5), 13.0; rtol=1e-12)
    @test isapprox(pm_limited_domain.mdot_flow_range.max(1.5), 18.5; rtol=1e-12)

    pm_bicubic = TM.TabulatedPerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
        [4.0 7.0 10.0; 6.0 9.0 12.0; 8.0 11.0 14.0],
        [0.70 0.80 0.90; 0.72 0.82 0.92; 0.74 0.84 0.94];
        interpolation=:bicubic,
    )
    vals_bicubic = TM.performance_from_stagnation(pm_bicubic, 2.5, 15.0, 300.0, 100_000.0)
    @test isapprox(vals_bicubic.pressure_ratio, 8.5; rtol=1e-12)

    from_inlet = TM.performance_from_stagnation(pm, 10_000.0, 15.0, 300.0, 100_000.0)
    @test isapprox(from_inlet.speed_coord, 10_000.0; rtol=1e-12)
    @test isapprox(from_inlet.flow_coord, 15.0; rtol=1e-12)
    @test isapprox(from_inlet.pressure_ratio, 4.5; rtol=1e-12)
    @test isapprox(from_inlet.eta, 0.91; rtol=1e-12)

    cmp_demo = TM.demo_tabulated_performance_map_compressor()
    cmp_vals = TM.performance_from_stagnation(
        cmp_demo,
        800.0,
        16.0,
        cmp_demo.Tt_ref,
        cmp_demo.Pt_ref,
    )
    @test cmp_vals.pressure_ratio > 1.0
    cmp_domain = TM.performance_map_domain(cmp_demo, cmp_demo.Tt_ref, cmp_demo.Pt_ref)
    @test cmp_domain.omega == (600.0, 1000.0)
    @test cmp_domain.mdot == (12.0, 20.0)
    @test isapprox(cmp_domain.mdot_flow_range.min(750.0), 12.0; rtol=1e-12)
    @test isapprox(cmp_domain.mdot_flow_range.max(750.0), 20.0; rtol=1e-12)

    @testset "Common Map Coordinate Helpers" begin
        @test isapprox(TM.corrected_speed(10_000.0, 300.0, 300.0), 10_000.0; rtol=1e-12)
        @test isapprox(TM.physical_speed(10_000.0, 300.0, 300.0), 10_000.0; rtol=1e-12)
        @test isapprox(TM.corrected_flow(15.0, 300.0, 100_000.0, 300.0, 100_000.0), 15.0; rtol=1e-12)
        @test isapprox(TM.physical_flow(15.0, 300.0, 100_000.0, 300.0, 100_000.0), 15.0; rtol=1e-12)
    end

    @testset "Linear Resampling" begin
        @test U.resample_linear([0.0, 1.0, 2.0], [0.0, 10.0, 20.0], 0.5) ≈ 5.0
        @test U.resample_linear([0.0, 1.0, 2.0], [0.0, 10.0, 20.0], [0.5, 1.5]) ≈ [5.0, 15.0]
        @test isnan(U.resample_linear([0.0, 1.0], [0.0, 10.0], -1.0; extrapolation=:nan))
        @test U.resample_linear([2.0, 1.0, 1.0, 0.0], [20.0, 10.0, 12.0, 0.0], 1.0; sort_inputs=true, deduplicate=:mean) ≈ 11.0
    end

    @testset "Turbomachine Residuals" begin
        eos = P.ideal_EOS()[:air]
        map = TM.TabulatedPerformanceMap(
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

        R_pout, R_dh_eff, R_P = TTM.compressor_residuals(
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

        R_pout, R_dh_eff, R_P = TTM.compressor_residuals(
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

        scales = TTM.compressor_residual_scales(
            pt_in,
            ht_in,
            pt_out,
            ht_out,
            mdot,
            omega,
            tau,
        )
        r_pout, r_dh_eff, r_P = TTM.compressor_residuals_scaled(
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

        r2_pout, r2_dh_eff, r2_P = TTM.compressor_residuals_scaled(
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

    @testset "Operating Point Solve (New API)" begin
        eos = P.ideal_EOS()[:air]
        demo_cmp = TP.Turbomachine.demo_tabulated_performance_map_compressor()
        demo_turb = TP.Turbomachine.demo_tabulated_performance_map_turbine()
        @test demo_cmp isa TP.Turbomachine.TabulatedPerformanceMap
        @test demo_turb isa TP.Turbomachine.TabulatedPerformanceMap

        common_vals = TP.Turbomachine.performance_from_stagnation(
            demo_cmp,
            800.0 * sqrt(432.0 / demo_cmp.Tt_ref),
            16.0 * (demo_cmp.Pt_ref / demo_cmp.Pt_ref) / sqrt(432.0 / demo_cmp.Tt_ref),
            432.0,
            demo_cmp.Pt_ref,
        )
        @test isapprox(common_vals.pressure_ratio, 1.80; rtol=1e-10)
        @test isapprox(common_vals.eta, 0.83; rtol=1e-10)
        @test isapprox(common_vals.speed_coord, 800.0; rtol=1e-10)
        @test isapprox(common_vals.flow_coord, 16.0; rtol=1e-10)
        common_domain = TP.Turbomachine.performance_map_domain(demo_cmp, 432.0, demo_cmp.Pt_ref)
        omega_scale = sqrt(432.0 / demo_cmp.Tt_ref)
        @test isapprox(common_domain.omega[1], 600.0 * omega_scale; rtol=1e-10)
        @test isapprox(common_domain.omega[2], 1000.0 * omega_scale; rtol=1e-10)

        map = TP.Turbomachine.TabulatedPerformanceMap(
            300.0,
            100_000.0,
            [1_000.0, 2_000.0],
            [10.0, 20.0],
            [1.8 2.2; 1.8 2.2],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
        )

        pt_in = 100_000.0
        Tt_in = 432.0
        ht_in = P.enthalpy_from_temperature(eos, Tt_in)
        pt_out = 200_000.0
        omega = 1_500.0 * sqrt(Tt_in / map.Tt_ref)

        mdot_expected = 15.0 * (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
        h2s = P.isentropic_enthalpy(eos, pt_in, ht_in, pt_out)
        ht_out_expected = ht_in + (h2s - ht_in) / 0.8
        tau_expected = mdot_expected * (ht_out_expected - ht_in) / omega

        sol = TP.Turbomachine.solve_operating_point(
            map,
            eos;
            pt_in=pt_in,
            ht_in=ht_in,
            pt_out=pt_out,
            omega=omega,
        )

        @test sol.converged
        @test sol.retcode == :success
        @test sol.n_candidates == 1
        candidate = only(sol.candidates)
        @test isapprox(candidate.mdot, mdot_expected; rtol=1e-8)
        @test isapprox(candidate.ht_out, ht_out_expected; rtol=1e-8)
        @test isapprox(candidate.tau, tau_expected; rtol=1e-8)
        @test isapprox(candidate.pressure_ratio, 2.0; rtol=1e-10)
        @test isapprox(candidate.efficiency, 0.8; rtol=1e-10)
        @test isapprox(candidate.residuals[1], 0.0; atol=1e-9)
        @test isapprox(candidate.residuals[2], 0.0; atol=1e-9)
        @test isapprox(candidate.residuals[3], 0.0; atol=1e-9)
        @test candidate.physically_admissible
        @test isapprox(candidate.machine_payload.speed_coord, 1_500.0; rtol=1e-10)
        @test isapprox(candidate.machine_payload.flow_coord, 15.0; rtol=1e-10)

        multi_root_map = TP.Turbomachine.TabulatedPerformanceMap(
            300.0,
            100_000.0,
            [1_000.0, 2_000.0],
            [10.0, 15.0, 20.0],
            [2.0 3.0 2.0; 2.0 3.0 2.0],
            [0.8 0.8 0.8; 0.8 0.8 0.8];
            interpolation=:bilinear,
        )

        multi_pt_in = 100_000.0
        multi_Tt_in = 432.0
        multi_ht_in = P.enthalpy_from_temperature(eos, multi_Tt_in)
        multi_pt_out = 250_000.0
        multi_omega = 1_500.0 * sqrt(multi_Tt_in / multi_root_map.Tt_ref)
        multi_flow_scale = (multi_pt_in / multi_root_map.Pt_ref) / sqrt(multi_Tt_in / multi_root_map.Tt_ref)

        multi_sol = TP.Turbomachine.solve_operating_point(
            multi_root_map,
            eos;
            pt_in=multi_pt_in,
            ht_in=multi_ht_in,
            pt_out=multi_pt_out,
            omega=multi_omega,
            continuation_hints=[12.0, 18.0],
        )

        @test multi_sol.converged
        @test multi_sol.retcode == :success
        @test multi_sol.n_candidates == 2
        @test isapprox(multi_sol.candidates[1].mdot, 12.5 * multi_flow_scale; atol=1e-8)
        @test isapprox(multi_sol.candidates[2].mdot, 17.5 * multi_flow_scale; atol=1e-8)

        low_branch = TP.Turbomachine.select_operating_point_branch(multi_sol; policy=:low)
        high_branch = TP.Turbomachine.select_operating_point_branch(multi_sol; policy=:high)
        near_low = TP.Turbomachine.select_operating_point_branch(
            multi_sol;
            policy=:nearest,
            reference_coordinate=12.0 * multi_flow_scale,
        )
        @test isapprox(low_branch.mdot, multi_sol.candidates[1].mdot; atol=1e-8)
        @test isapprox(high_branch.mdot, multi_sol.candidates[2].mdot; atol=1e-8)
        @test isapprox(near_low.mdot, multi_sol.candidates[1].mdot; atol=1e-8)

        no_solution = TP.Turbomachine.solve_operating_point(
            map,
            eos;
            pt_in=pt_in,
            ht_in=ht_in,
            pt_out=350_000.0,
            omega=omega,
        )
        @test !no_solution.converged
        @test no_solution.retcode == :no_candidate
        @test isempty(no_solution.candidates)

        axial_model = TP.Turbomachine.demo_axial_compressor_model()
        axial_map = TP.Turbomachine.tabulate_axial_machine_model(
            axial_model;
            n_speed=5,
            n_flow=7,
            Tt_in_ref=288.15,
            Pt_in_ref=101_325.0,
            Tt_ref=288.15,
            Pt_ref=101_325.0,
            interpolation=:bilinear,
            boundary_resolution=81,
            want_diagnostics=false,
        )
        omega_diag = 800.0
        mdot_diag = 16.0
        axial_vals = TP.Turbomachine.performance_from_stagnation(
            axial_map,
            omega_diag,
            mdot_diag,
            axial_map.Tt_ref,
            axial_map.Pt_ref,
        )
        diag_sol = TP.Turbomachine.solve_operating_point(
            axial_map,
            eos;
            pt_in=axial_map.Pt_ref,
            ht_in=P.enthalpy_from_temperature(eos, axial_map.Tt_ref),
            pt_out=axial_vals.pressure_ratio * axial_map.Pt_ref,
            omega=omega_diag,
        )
        @test diag_sol.converged
        @test diag_sol.n_candidates >= 1
        diag_candidate = TP.Turbomachine.select_operating_point_branch(
            diag_sol;
            policy=:nearest,
            reference_coordinate=mdot_diag,
        )
        @test diag_candidate !== nothing
        replay = TP.Turbomachine.replay_operating_point_with_streamtube(
            axial_model,
            eos,
            diag_candidate;
            pt_in=axial_map.Pt_ref,
            ht_in=P.enthalpy_from_temperature(eos, axial_map.Tt_ref),
            omega=omega_diag,
        )
        @test length(replay.physical.station_data) == length(axial_model.rows) + 1
        @test length(replay.physical.row_data) == length(axial_model.rows)
        @test haskey(replay.physical.station_data[1], :pt_t)
        @test hasproperty(replay.physical.row_data[1], :aero)
        @test hasproperty(replay.physical.row_data[1], :choke)
        @test hasproperty(replay.physical.row_data[1], :outlet_candidates)
        @test hasproperty(replay.physical.row_data[1], :selected_candidate_index)
        @test hasproperty(replay.physical.row_data[1].aero, :incidence)
        @test hasproperty(replay.physical.row_data[1].aero, :delta_s_t)
        @test haskey(replay.physical.station_data[1], :mdot_station)
        @test haskey(replay.summary, :pressure_ratio_error)
        @test haskey(replay.summary, :ht_out_error)
        @test haskey(replay.summary, :streamtube_thermo_efficiency)

        direct_diag = TP.Turbomachine.diagnose_axial_operating_point(
            axial_model,
            eos;
            pt_in=axial_map.Pt_ref,
            ht_in=P.enthalpy_from_temperature(eos, axial_map.Tt_ref),
            omega=omega_diag,
            Vx_inlet=inlet_Vx(
                axial_model,
                eos,
                axial_map.Pt_ref,
                P.enthalpy_from_temperature(eos, axial_map.Tt_ref),
                mdot_diag,
                0.0,
            ),
        )
        @test haskey(direct_diag.summary, :pressure_ratio)
        @test haskey(direct_diag.summary, :power)
        @test haskey(direct_diag.summary, :thermo_efficiency)
        @test length(direct_diag.physical.row_data) == length(axial_model.rows)

        turbine_model = TP.Turbomachine.demo_axial_turbine_model()
        turbine_diag = TP.Turbomachine.diagnose_axial_operating_point(
            turbine_model,
            eos;
            pt_in=101_325.0,
            ht_in=P.enthalpy_from_temperature(eos, 288.15),
            omega=600.0,
            Vx_inlet=inlet_Vx(
                turbine_model,
                eos,
                101_325.0,
                P.enthalpy_from_temperature(eos, 288.15),
                8.0,
                0.0,
            ),
        )
        @test turbine_diag.summary.pressure_ratio < 1.0
        @test isapprox(turbine_diag.summary.efficiency, turbine_diag.summary.thermo_efficiency; rtol=1e-10, atol=1e-10)
    end

    @testset "Operating Point Sweep (New API)" begin
        eos = P.ideal_EOS()[:air]
        map = TP.Turbomachine.TabulatedPerformanceMap(
            300.0,
            100_000.0,
            [1_000.0, 2_000.0],
            [10.0, 20.0],
            [1.8 2.2; 1.8 2.2],
            [0.8 0.8; 0.8 0.8];
            interpolation=:bilinear,
        )

        pt_in = 100_000.0
        Tt_in = 432.0
        ht_in = P.enthalpy_from_temperature(eos, Tt_in)
        omega_grid = sqrt(Tt_in / map.Tt_ref) .* [1_000.0, 1_500.0, 2_000.0]
        pt_out_grid = [200_000.0]
        sweep = TP.Turbomachine.solve_operating_sweep(
            map,
            eos;
            pt_in=pt_in,
            ht_in=ht_in,
            omega_grid=omega_grid,
            pt_out_grid=pt_out_grid,
            branch_selection=:nearest,
            initial_branch_coordinate=15.0,
        )

        @test sweep.mode == :single
        @test size(sweep.mdot) == (length(omega_grid), length(pt_out_grid))
        @test all(sweep.converged)
        for i in eachindex(omega_grid)
            mdot_expected = 15.0 * (pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
            @test isapprox(sweep.mdot[i, 1], mdot_expected; rtol=1e-8)
            @test isapprox(sweep.PR[i, 1], 2.0; rtol=1e-10)
            @test isapprox(sweep.eta[i, 1], 0.8; rtol=1e-10)
        end

        multi_root_map = TP.Turbomachine.TabulatedPerformanceMap(
            300.0,
            100_000.0,
            [1_000.0, 2_000.0],
            [10.0, 15.0, 20.0],
            [2.0 3.0 2.0; 2.0 3.0 2.0],
            [0.8 0.8 0.8; 0.8 0.8 0.8];
            interpolation=:bilinear,
        )
        sweep_all = TP.Turbomachine.solve_operating_sweep(
            multi_root_map,
            eos;
            pt_in=pt_in,
            ht_in=ht_in,
            omega_grid=fill(1_500.0 * sqrt(Tt_in / multi_root_map.Tt_ref), 2),
            pt_out_grid=[250_000.0],
            track_all_branches=true,
        )

        @test sweep_all.mode == :all
        @test length(sweep_all.rows) == 4
        good_rows = filter(row -> row.converged, sweep_all.rows)
        @test length(good_rows) == 4
        @test length(unique(row.branch_id for row in good_rows)) == 2
        @test all(row -> isfinite(row.mdot), good_rows)
    end

    @testset "Map IO" begin
        common_map = TP.Turbomachine.demo_tabulated_performance_map_compressor(; interpolation=:bilinear)
        common_map_toml_path = tempname() * ".toml"
        TP.Turbomachine.write_toml(common_map, common_map_toml_path)
        common_map_loaded = TP.Turbomachine.read_toml(TP.Turbomachine.TabulatedPerformanceMap, common_map_toml_path)
        @test common_map_loaded isa TP.Turbomachine.TabulatedPerformanceMap
        @test TP.Turbomachine.tabulated_speed_grid(common_map_loaded) == TP.Turbomachine.tabulated_speed_grid(common_map)
        @test TP.Turbomachine.tabulated_flow_grid(common_map_loaded) == TP.Turbomachine.tabulated_flow_grid(common_map)
        @test TP.Turbomachine.tabulated_pressure_ratio_table(common_map_loaded) == TP.Turbomachine.tabulated_pressure_ratio_table(common_map)
        @test TP.Turbomachine.tabulated_eta_table(common_map_loaded) == TP.Turbomachine.tabulated_eta_table(common_map)

        table_map_toml_path = tempname() * ".toml"
        U.write_toml(pm.pressure_ratio_map, table_map_toml_path)
        table_map_toml_loaded = U.read_toml(U.AbstractTableMap, table_map_toml_path)
        @test U.table_interpolation(table_map_toml_loaded) == :bilinear
        @test U.table_xgrid(table_map_toml_loaded) == U.table_xgrid(pm.pressure_ratio_map)
        @test U.table_ygrid(table_map_toml_loaded) == U.table_ygrid(pm.pressure_ratio_map)
        @test U.table_values(table_map_toml_loaded) == U.table_values(pm.pressure_ratio_map)

        vals_ref = TM.performance_from_stagnation(pm_bicubic, 2.5, 15.0, 300.0, 100_000.0)
        compressor_map_toml_path = tempname() * ".toml"
        TM.write_toml(pm_bicubic, compressor_map_toml_path)
        pm_toml_loaded = TM.read_toml(TM.TabulatedPerformanceMap, compressor_map_toml_path)
        vals_toml_loaded = TM.performance_from_stagnation(pm_toml_loaded, 2.5, 15.0, 300.0, 100_000.0)
        @test isapprox(vals_toml_loaded.pressure_ratio, vals_ref.pressure_ratio; rtol=1e-12)
        @test isapprox(vals_toml_loaded.eta, vals_ref.eta; rtol=1e-12)

        meanline = TTM.demo_axial_compressor_model()
        m_mid = 0.5 * (meanline.m_tip_bounds[1] + meanline.m_tip_bounds[2])
        Vx_mid = 0.5 * (meanline.Vx_bounds[1] + meanline.Vx_bounds[2])
        radii = AM.meanline_radii(meanline)
        eos_axial = P.IdealGasEOS(:axial; gas_constant=meanline.gas_constant, gamma=meanline.gamma)
        Tt_ref = 300.0
        Pt_ref = 101_325.0
        ht_ref = P.enthalpy_from_temperature(eos_axial, Tt_ref)
        a0 = P.speed_of_sound(eos_axial, Pt_ref, ht_ref)
        omega_mid = m_mid * a0 / meanline.r_tip_ref
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
        meanline_vals = AM.streamtube_solve(
            meanline,
            eos_axial,
            radii,
            Pt_ref,
            ht_ref,
            omega_mid,
            Vx_mid,
            0.0,
        )
        @test isfinite(meanline_vals.PR)
        @test isfinite(meanline_vals.eta)
        meanline_vals_core = AM.streamtube_solve(
            meanline,
            eos_axial,
            radii,
            Pt_ref,
            ht_ref,
            omega_mid,
            Vx_mid,
            0.0,
        )
        @test isfinite(meanline_vals_core.PR)
        @test isfinite(meanline_vals_core.eta)
        @test isapprox(meanline_vals_core.PR, meanline_vals.PR; rtol=1e-12)
        @test isapprox(meanline_vals_core.eta, meanline_vals.eta; rtol=1e-12)
        meanline_vals_from_mdot = AM.streamtube_solve(
            meanline,
            eos_axial,
            radii,
            Pt_ref,
            ht_ref,
            omega_mid,
            inlet_Vx(meanline, eos_axial, Pt_ref, ht_ref, meanline_vals_core.stations[1].mdot_station, 0.0),
            0.0,
        )
        @test isapprox(meanline_vals_from_mdot.PR, meanline_vals_core.PR; rtol=1e-10)
        @test isapprox(meanline_vals_from_mdot.eta, meanline_vals_core.eta; rtol=1e-9)
        @test isapprox(meanline_vals_from_mdot.stations[1].Vx, meanline_vals_core.stations[1].Vx; rtol=1e-9)
        feasible = AM.feasible_flow_limits(
            meanline,
            eos_axial,
            [omega_mid],
            meanline.Vx_bounds[1],
            meanline.Vx_bounds[2];
            pt_in=Pt_ref,
            ht_in=ht_ref,
        )
        @test feasible.valid_speed_idx == [1]
        @test length(feasible.Vx_min) == 1
        @test feasible.Vx_max[1] >= feasible.Vx_min[1]
        @test isfinite(feasible.mdot_min[1])
        @test isfinite(feasible.mdot_max[1])
        @test feasible.mdot_max[1] >= feasible.mdot_min[1]
        sampled = AM.sample_streamtube_solve(
            meanline,
            eos_axial,
            [omega_mid],
            [Vx_mid];
            pt_in=Pt_ref,
            ht_in=ht_ref,
            Vx_min=feasible.Vx_min,
            Vx_max=feasible.Vx_max,
        )
        @test size(sampled.pr_table) == (1, 1)
        @test size(sampled.eta_table) == (1, 1)
        @test size(sampled.mdot_table) == (1, 1)
        @test isfinite(sampled.mdot_table[1, 1])
        @test length(meanline_vals_core.station_data) == length(meanline.rows) + 1
        @test length(meanline_vals_core.row_data) == length(meanline.rows)
        @test meanline_vals_core.station_data[1].area ≈ AM.station_area(meanline, 1)
        @test AM.row_annulus_area(meanline.rows[1]) > 0
        @test AM.row_mean_radius(meanline.rows[1]) > meanline.rows[1].r_hub
        @test AM.row_angular_speed(meanline.rows[1], omega_mid) == meanline.rows[1].speed_ratio_to_ref * omega_mid

        dim_from_meanline = TTM.tabulate_compressor_meanline_model_dim(
            meanline;
            n_speed=9,
            n_flow=11,
            interpolation=:bilinear,
            Tt_in_ref=300.0,
            Pt_in_ref=101_325.0,
            Tt_ref=300.0,
            Pt_ref=101_325.0,
            eos=eos_axial,
        )
        @test dim_from_meanline isa TTM.TabulatedPerformanceMap
        dim_grid_omega = collect(TTM.tabulated_speed_grid(dim_from_meanline))
        dim_grid_mdot_corr = collect(TTM.tabulated_flow_grid(dim_from_meanline))
        @test length(dim_grid_mdot_corr) == 11
        @test length(dim_grid_omega) >= 2
        i_mid = cld(length(dim_grid_omega), 2)
        j_mid = cld(length(dim_grid_mdot_corr), 2)
        vals_dim = TTM.performance_from_stagnation(
            dim_from_meanline,
            TTM.physical_speed(dim_grid_omega[i_mid], 300.0, 300.0),
            TTM.physical_flow(dim_grid_mdot_corr[j_mid], 300.0, 101_325.0, 300.0, 101_325.0),
            300.0,
            101_325.0,
        )
        omega_sample = TTM.physical_speed(dim_grid_omega[i_mid], 300.0, 300.0)
        mdot_sample = TTM.physical_flow(dim_grid_mdot_corr[j_mid], 300.0, 101_325.0, 300.0, 101_325.0)
        Vx_sample = inlet_Vx(meanline, eos_axial, 101_325.0, ht_ref, mdot_sample, 0.0)
        vals_streamtube = AM.streamtube_solve(
            meanline,
            eos_axial,
            radii,
            101_325.0,
            ht_ref,
            omega_sample,
            Vx_sample,
            0.0,
        )
        @test isapprox(vals_dim.pressure_ratio, vals_streamtube.PR; rtol=1e-8)
        @test isapprox(vals_dim.eta, vals_streamtube.eta; rtol=1e-8)

        shared_from_meanline = TP.Turbomachine.tabulate_axial_machine_model(
            meanline;
            eos=eos_axial,
            n_speed=7,
            n_flow=9,
            interpolation=:bilinear,
            Tt_in_ref=300.0,
            Pt_in_ref=101_325.0,
            Tt_ref=300.0,
            Pt_ref=101_325.0,
            boundary_resolution=101,
        )
        @test shared_from_meanline isa TP.Turbomachine.TabulatedPerformanceMap
        shared_omega_grid = TP.Turbomachine.tabulated_speed_grid(shared_from_meanline)
        shared_mdot_grid = TP.Turbomachine.tabulated_flow_grid(shared_from_meanline)
        @test length(shared_omega_grid) >= 2
        @test length(shared_mdot_grid) == 9
        vals_shared = TP.Turbomachine.performance_from_stagnation(
            shared_from_meanline,
            shared_omega_grid[cld(length(shared_omega_grid), 2)],
            shared_mdot_grid[cld(length(shared_mdot_grid), 2)],
            300.0,
            101_325.0,
        )
        @test isfinite(vals_shared.pressure_ratio)
        @test isfinite(vals_shared.eta)
        shared_domain = TP.Turbomachine.performance_map_domain(
            shared_from_meanline,
            300.0,
            101_325.0,
        )
        @test shared_domain.omega[1] > 0
        @test shared_domain.omega[2] > shared_domain.omega[1]

        meanline_path = tempname() * ".toml"
        AM.write_toml(meanline, meanline_path)
        meanline_loaded = AM.read_toml(TTM.AxialModel, meanline_path)
        @test_throws ErrorException TM.read_toml(TM.TabulatedPerformanceMap, meanline_path; group="axial_model")
        Vx_mid_loaded = 0.5 * (meanline_loaded.Vx_bounds[1] + meanline_loaded.Vx_bounds[2])
        meanline_vals_loaded = AM.streamtube_solve(
            meanline_loaded,
            eos_axial,
            radii,
            Pt_ref,
            ht_ref,
            omega_mid,
            Vx_mid_loaded,
            0.0,
        )
        @test isapprox(meanline_vals_loaded.PR, meanline_vals.PR; rtol=1e-12)
        @test isapprox(meanline_vals_loaded.eta, meanline_vals.eta; rtol=1e-12)

        meanline_turbine = TTM.demo_axial_turbine_model()
        m_mid_t = 0.5 * (meanline_turbine.m_tip_bounds[1] + meanline_turbine.m_tip_bounds[2])
        Vx_mid_t = 0.5 * (meanline_turbine.Vx_bounds[1] + meanline_turbine.Vx_bounds[2])
        eos_turbine = P.IdealGasEOS(:axial; gas_constant=meanline_turbine.gas_constant, gamma=meanline_turbine.gamma)
        ht_ref_t = P.enthalpy_from_temperature(eos_turbine, Tt_ref)
        a0_t = P.speed_of_sound(eos_turbine, Pt_ref, ht_ref_t)
        omega_mid_t = m_mid_t * a0_t / meanline_turbine.r_tip_ref
        vals_turbine = AM.streamtube_solve(
            meanline_turbine,
            eos_turbine,
            AM.meanline_radii(meanline_turbine),
            Pt_ref,
            ht_ref_t,
            omega_mid_t,
            Vx_mid_t,
            0.0,
        )
        @test hasproperty(vals_turbine, :valid)
        @test hasproperty(vals_turbine, :PR)
        @test hasproperty(vals_turbine, :eta)

        turbine_map = TM.demo_tabulated_performance_map_turbine(; interpolation=:bilinear)
        vals_t_ref = TM.performance_from_stagnation(
            turbine_map,
            TM.physical_speed(0.8, turbine_map.Tt_ref, turbine_map.Tt_ref),
            TM.physical_flow(1.8, turbine_map.Tt_ref, turbine_map.Pt_ref, turbine_map.Tt_ref, turbine_map.Pt_ref),
            turbine_map.Tt_ref,
            turbine_map.Pt_ref,
        )
        turbine_map_toml_path = tempname() * ".toml"
        TM.write_toml(turbine_map, turbine_map_toml_path)
        turbine_map_toml_loaded = TM.read_toml(TM.TabulatedPerformanceMap, turbine_map_toml_path)
        vals_t_toml_loaded = TM.performance_from_stagnation(
            turbine_map_toml_loaded,
            TM.physical_speed(0.8, turbine_map.Tt_ref, turbine_map.Tt_ref),
            TM.physical_flow(1.8, turbine_map.Tt_ref, turbine_map.Pt_ref, turbine_map.Tt_ref, turbine_map.Pt_ref),
            turbine_map.Tt_ref,
            turbine_map.Pt_ref,
        )
        @test isapprox(vals_t_toml_loaded.pressure_ratio, vals_t_ref.pressure_ratio; rtol=1e-12)
        @test isapprox(vals_t_toml_loaded.eta, vals_t_ref.eta; rtol=1e-12)
    end
end
