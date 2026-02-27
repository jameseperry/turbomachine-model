module Utility

export AbstractTableMap
export BilinearMap, BicubicMap
export interpolation_map
export bilinear_evaluate, bicubic_evaluate, table_evaluate
export table_xgrid, table_ygrid, table_values, table_interpolation
export write_toml, read_toml

include("bilinear_map.jl")
include("bicubic_map.jl")
include("toml_io.jl")

table_xgrid(map::AbstractTableMap) = error("table_xgrid not implemented for $(typeof(map))")
table_ygrid(map::AbstractTableMap) = error("table_ygrid not implemented for $(typeof(map))")
table_values(map::AbstractTableMap) = error("table_values not implemented for $(typeof(map))")
table_interpolation(map::AbstractTableMap) = error("table_interpolation not implemented for $(typeof(map))")

function interpolation_map(
    interpolation::Symbol,
    xgrid::Vector{<:Real},
    ygrid::Vector{<:Real},
    table::Matrix{<:Real},
)
    return interpolation_map(Val(interpolation), xgrid, ygrid, table)
end

interpolation_map(
    ::Val{:bilinear},
    xgrid::Vector{<:Real},
    ygrid::Vector{<:Real},
    table::Matrix{<:Real},
) = BilinearMap(xgrid, ygrid, table)

interpolation_map(
    ::Val{:bicubic},
    xgrid::Vector{<:Real},
    ygrid::Vector{<:Real},
    table::Matrix{<:Real},
) = BicubicMap(xgrid, ygrid, table)

function interpolation_map(
    ::Val{K},
    xgrid::Vector{<:Real},
    ygrid::Vector{<:Real},
    table::Matrix{<:Real},
) where {K}
    error("unsupported interpolation=$(K) (expected :bilinear or :bicubic)")
end

end # module Utility
