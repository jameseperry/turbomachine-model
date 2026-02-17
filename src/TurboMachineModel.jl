module TurboMachineModel

include("framework/Framework.jl")
include("physics/Physics.jl")
include("components/Components.jl")

using .Framework
using .Physics
using .Components

export Framework, Physics, Components

end # module TurboMachineModel
