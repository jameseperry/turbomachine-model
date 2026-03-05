module Turbomachine

abstract type AbstractCompressorPerformanceMap end

include("nondimensional_performance_map.jl")
include("axial_machine/AxialMachine.jl")
include("compressor/Compressor.jl")
include("turbine/Turbine.jl")

import .AxialMachine
import .Compressor
import .Turbine

export NondimensionalPerformanceMap
export AxialMachine
export Compressor
export Turbine

end # module Turbomachine
