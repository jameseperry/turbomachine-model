module Physics

include("fluids/Fluids.jl")
include("turbomachine/Turbomachine.jl")

using .Fluids
using .Turbomachine

export Fluids
export Turbomachine

end # module Physics
