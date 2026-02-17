module Mechanical

using ..AbstractComponentDef: AbstractComponent
import ..Components: port_specs, required_ports, validate
using ..PortPresets: SHAFT_PORT

include("inertial_shaft.jl")

export InertialShaft

"""
Mechanical-side component namespace.
Component data types are intentionally minimal at this stage.
"""

end # module Mechanical
