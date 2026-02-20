@testset "Fluid Properties" begin
    P = TurboMachineModel.Physics.Fluids

    @testset "Ideal Gas Composition" begin
        comp = P.IdealGasComposition(:air; gas_constant=287.05, gamma=1.4)
        @test comp.id == :air
        @test isapprox(comp.gamma, 1.4; rtol=1e-12)
        @test isapprox(comp.specific_heat_cp - comp.specific_heat_cv, comp.gas_constant; rtol=1e-12)

        T = 300.0
        h = P.enthalpy_from_temperature(comp, T)
        @test isapprox(P.temperature_from_enthalpy(comp, h), T; rtol=1e-12)

        p = 101_325.0
        rho = P.density_from_pressure_temperature(comp, p, T)
        @test isapprox(P.density(comp, p, h), rho; rtol=1e-12)

        a = P.speed_of_sound_from_temperature(comp, T)
        @test isapprox(P.speed_of_sound(comp, p, h), a; rtol=1e-12)
    end

    @testset "Equations of State Registry" begin
        laws = P.ideal_EOS()
        @test :air in keys(laws)
        @test :steam in keys(laws)

        p = 101_325.0
        h_air = 300_000.0
        T_ph = P.temperature(laws, :air, p, h_air)
        s_ph = P.entropy(laws, :air, p, h_air)
        @test T_ph > 0
        @test P.density(laws, :air, p, h_air) > 0
        @test P.speed_of_sound(laws, :air, p, h_air) > 0
        @test P.heat_capacity_cp(laws, :air, p, h_air) > 0
        @test P.dynamic_viscosity(laws, :air, p, h_air) > 0
        @test P.thermal_conductivity(laws, :air, p, h_air) > 0
        @test P.phase(laws, :air, p, h_air) == :gas
        @test isapprox(P.enthalpy_from_pressure_entropy(laws, :air, p, s_ph), h_air; rtol=1e-12)

        @test_throws ErrorException P.temperature(laws, :unknown, p, h_air)
    end
end
