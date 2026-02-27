using HDF5

function _write_hdf5_group!(group::Union{HDF5.File, HDF5.Group}, map::AbstractTableMap)
    write(group, "xgrid", table_xgrid(map))
    write(group, "ygrid", table_ygrid(map))
    write(group, "table", table_values(map))
    attrs(group)["interpolation"] = String(table_interpolation(map))
    return nothing
end

function _read_hdf5_group(group::Union{HDF5.File, HDF5.Group})
    haskey(group, "xgrid") || error("missing dataset xgrid")
    haskey(group, "ygrid") || error("missing dataset ygrid")
    haskey(group, "table") || error("missing dataset table")
    haskey(attrs(group), "interpolation") || error("missing attribute interpolation")

    interpolation_raw = attrs(group)["interpolation"]
    interpolation = Symbol(String(interpolation_raw))
    xgrid = vec(Float64.(read(group["xgrid"])))
    ygrid = vec(Float64.(read(group["ygrid"])))
    table = Float64.(read(group["table"]))
    return interpolation_map(interpolation, xgrid, ygrid, table)
end

function write_hdf5(
    parent::Union{HDF5.File, HDF5.Group},
    name::AbstractString,
    map::AbstractTableMap,
)
    haskey(parent, name) && error("HDF5 object $(name) already exists")
    group = create_group(parent, name)
    _write_hdf5_group!(group, map)
    return nothing
end

function read_hdf5(
    ::Type{AbstractTableMap},
    parent::Union{HDF5.File, HDF5.Group},
    name::AbstractString,
)
    haskey(parent, name) || error("missing HDF5 object $(name)")
    return _read_hdf5_group(parent[name])
end

function write_hdf5(
    map::AbstractTableMap,
    path::AbstractString;
    group::AbstractString="table_map",
)
    h5open(path, "w") do file
        write_hdf5(file, group, map)
    end
    return path
end

function read_hdf5(
    ::Type{AbstractTableMap},
    path::AbstractString;
    group::AbstractString="table_map",
)
    return h5open(path, "r") do file
        read_hdf5(AbstractTableMap, file, group)
    end
end
