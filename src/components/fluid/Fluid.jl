module Fluid

using ..AbstractComponentDef: AbstractComponent
import ..Components: port_specs, required_ports, validate
using ..PortPresets: FLUID_PORT, SHAFT_PORT

include("thermo_utils.jl")
include("performance_map.jl")
include("combustor.jl")
include("turbo_machine_section.jl")

export TurboMachineSection, Combustor
export static_temperature_from_total, total_temperature_from_static
export static_pressure_from_total, total_pressure_from_static
export static_enthalpy_from_total, total_enthalpy_from_static
export PerformanceMap, corrected_speed, corrected_flow
export map_pr_eta, map_pr_eta_from_inlet

"""
Fluid-side component namespace.
Component data types are intentionally minimal at this stage.
"""

end # module Fluid
