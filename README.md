# TurboMachineModel

Julia project scaffold for turbomachinery simulation.

Scope for future implementation:
- Gas-side modeling (thermodynamics/flow)
- Mechanical-side modeling (rotordynamics/components)
- Coupling/integration workflows

## Getting started

1. Install Julia 1.10+.
2. From this directory, start Julia and activate the environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

3. Run tests:

```julia
using Pkg
Pkg.test()
```
