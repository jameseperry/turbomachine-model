module Physics

include("fluids/Fluids.jl")
include("turbomachine.jl")

using .Fluids

export Fluids
export turbomachine_residuals

end # module Physics
