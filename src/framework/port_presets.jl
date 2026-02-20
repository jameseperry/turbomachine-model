const FLUID_THROUGH_VARS = [
    (var=:pt, unit=:Pa, coupling=:equality),
    (var=:ht, unit=:J_per_kg, coupling=:equality),
    (var=:mdot, unit=:kg_per_s, coupling=:conservation),
    (var=:composition, unit=:symbol, coupling=:equality),
]

const SHAFT_VARS = [
    (var=:omega, unit=:rad_per_s, coupling=:equality),
    (var=:tau, unit=:N_m, coupling=:conservation),
]

const THERMAL_VARS = [
    (var=:temperature, unit=:K, coupling=:equality),
    (var=:heat_flow, unit=:W, coupling=:conservation),
]

const FLUID_PORT = (domain=:fluid, vars=FLUID_THROUGH_VARS)
const SHAFT_PORT = (domain=:shaft, vars=SHAFT_VARS)
const THERMAL_PORT = (domain=:thermal, vars=THERMAL_VARS)
