module AxialMachine

export AbstractRowAeroModel
export RotorAeroModel, StatorAeroModel
export RowAeroOutput
export AxialRow, AxialMachineModel
export row_aero
export streamtube_solve
export sample_streamtube_solve
export feasible_flow_limits

include("row_aero_model.jl")
include("axial_machine_model.jl")
include("axial_machine_model_io.jl")
include("streamtube_solver.jl")
include("flow_limits.jl")

end # module AxialMachine
