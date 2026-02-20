const FLUID_THROUGH_VARS = [
    (var=:p, unit=:Pa, coupling=:equality),
    (var=:h, unit=:J_per_kg, coupling=:equality),
    (var=:mdot, unit=:kg_per_s, coupling=:conservation),
    (var=:fluid_id, unit=:dimensionless, coupling=:equality),
]

const SHAFT_VARS = [
    (var=:omega, unit=:rad_per_s, coupling=:equality),
    (var=:tau, unit=:N_m, coupling=:conservation),
]

const FLUID_PORT = (domain=:fluid, vars=FLUID_THROUGH_VARS)
const SHAFT_PORT = (domain=:shaft, vars=SHAFT_VARS)
