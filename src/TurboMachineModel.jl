module TurboMachineModel

include("components/Components.jl")
include("structure/Structure.jl")

using .Components
using .Structure

export Structure, Components

end # module TurboMachineModel
