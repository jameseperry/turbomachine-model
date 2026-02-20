using Test
using TurboMachineModel

@testset "TurboMachineModel.jl" begin
    include("framework_smoke_tests.jl")
    include("component_tests.jl")
    include("fluid_properties_tests.jl")
    include("physics_tests.jl")
    include("network_tests.jl")
end
