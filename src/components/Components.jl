module Components

import ..Framework: AbstractComponent, port_specs, required_ports, validate
using ..Framework: FLUID_PORT, SHAFT_PORT
using ..Physics.Fluids: PerformanceMap

include("combustor.jl")
include("plenum.jl")
include("turbo_machine_section.jl")
include("inertial_shaft.jl")
include("gearbox.jl")

export AbstractComponent, Combustor, Plenum, TurboMachineSection, InertialShaft, Gearbox

"""
Namespace for concrete component implementations
(compressor/turbine sections, combustors, shafts, etc.).
"""

end # module Components
