module TurboMachineModel

include("structure/Structure.jl")
include("components/Components.jl")

using .Structure
using .Components

export Structure, Components

end # module TurboMachineModel
