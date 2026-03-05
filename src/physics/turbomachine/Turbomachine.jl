module Turbomachine

include("axial_machine/AxialMachine.jl")
include("compressor/Compressor.jl")
include("turbine/Turbine.jl")

import .AxialMachine
import .Compressor
import .Turbine

export AxialMachine
export Compressor
export Turbine

end # module Turbomachine
