@testset "Framework Smoke" begin
    @test isdefined(TurboMachineModel, :Framework)
    @test isdefined(TurboMachineModel, :Physics)
    @test isdefined(TurboMachineModel, :Components)
    @test isdefined(TurboMachineModel.Framework, :Network)
    @test isdefined(TurboMachineModel.Framework, :PortSpec)
    @test isdefined(TurboMachineModel.Components, :AbstractComponent)
    @test isdefined(TurboMachineModel.Components, :Combustor)
    @test isdefined(TurboMachineModel.Components, :TurboMachineSection)
    @test isdefined(TurboMachineModel.Components, :InertialShaft)
end
