using Test
using TurboMachineModel

@testset "TurboMachineModel.jl" begin
    @test isdefined(TurboMachineModel, :Structure)
    @test isdefined(TurboMachineModel, :Components)
end
