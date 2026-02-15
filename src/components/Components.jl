module Components

include("fluid/Fluid.jl")
include("mechanical/Mechanical.jl")

using .Fluid
using .Mechanical

export Fluid, Mechanical

"""
Namespace for concrete component implementations
(compressor, turbine, combustor, shaft, etc.).
Implementation intentionally deferred.
"""

end # module Components
