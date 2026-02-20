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
        @test isapprox(P.density_from_pressure_enthalpy(comp, p, h), rho; rtol=1e-12)

        a = P.speed_of_sound_from_temperature(comp, T)
        @test isapprox(P.speed_of_sound_from_enthalpy(comp, h), a; rtol=1e-12)
    end

    @testset "Equations of State Registry" begin
        laws = P.ideal_EOS()
        @test :air in keys(laws)
        @test :steam in keys(laws)

        p = 101_325.0
        T = 300.0
        h_air = P.enthalpy_from_temperature(laws, :air, T)
        @test isapprox(P.temperature_from_enthalpy(laws, :air, h_air), T; rtol=1e-12)
        @test P.density_from_pressure_temperature(laws, :air, p, T) > 0
        @test P.speed_of_sound_from_temperature(laws, :air, T) > 0

        @test_throws ErrorException P.temperature_from_enthalpy(laws, :unknown, h_air)
    end
end
