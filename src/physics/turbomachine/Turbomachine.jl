module Turbomachine

include("compressor/Compressor.jl")
include("turbine/Turbine.jl")

using .Compressor
using .Turbine

export Compressor
export Turbine

end # module Turbomachine
