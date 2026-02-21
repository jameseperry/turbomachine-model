module Components

import ..Component: AbstractComponent, ports, variables
using ..Component: ComponentPort, ComponentVariable
using ..Physics.Fluids: PerformanceMap

include("combustor.jl")
include("plenum.jl")
include("turbo_machine_section.jl")
include("inertial_shaft.jl")
include("gearbox.jl")

export Combustor, Plenum, TurboMachineSection, InertialShaft, Gearbox

"""
Namespace for concrete component implementations
(compressor/turbine sections, combustors, shafts, etc.).
"""

end # module Components
