@testset "Fluid Properties" begin
    P = TurboMachineModel.Physics.Fluids

    @testset "Ideal Gas Composition" begin
        comp = P.IdealGasEOS(:air; gas_constant=287.05, gamma=1.4)
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

        p2 = 2 * p
        h2s = P.isentropic_enthalpy(comp, p, h, p2)
        s1 = P.entropy(comp, p, h)
        s2 = P.entropy(comp, p2, h2s)
        @test isapprox(s2, s1; rtol=1e-12)
    end

    @testset "Equations of State Registry" begin
        laws = P.ideal_EOS()
        @test :air in keys(laws)
        @test :steam in keys(laws)
        @test laws[:air] isa P.AbstractEOS
        @test laws[:steam] isa P.AbstractEOS
        @test_throws KeyError laws[:unknown]

        air = laws[:air]

        p = 101_325.0
        h_air = 300_000.0
        T_ph = P.temperature(air, p, h_air)
        s_ph = P.entropy(air, p, h_air)
        @test T_ph > 0
        @test P.density(air, p, h_air) > 0
        @test P.speed_of_sound(air, p, h_air) > 0
        @test P.heat_capacity_cp(air, p, h_air) > 0
        @test P.dynamic_viscosity(air, p, h_air) > 0
        @test P.thermal_conductivity(air, p, h_air) > 0
        @test P.phase(air, p, h_air) == :gas
        @test isapprox(P.enthalpy_from_pressure_entropy(air, p, s_ph), h_air; rtol=1e-12)
        h2s = P.isentropic_enthalpy(air, p, h_air, 2 * p)
        s2 = P.entropy(air, 2 * p, h2s)
        @test isapprox(s2, s_ph; rtol=1e-12)
    end
end
