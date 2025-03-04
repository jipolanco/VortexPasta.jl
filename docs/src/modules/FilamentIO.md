# FilamentIO

```@meta
CurrentModule = VortexPasta.FilamentIO
CollapsedDocStrings = true
```

```@docs
FilamentIO
```

## VTKHDF format

The [VTKHDF format](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format)
can be used to store data from numerical simulations and visualise it using ParaView.
It is based on HDF5, meaning that data can be accessed using the HDF5 libraries and tools.

### Writing data

```@docs
write_vtkhdf
Base.setindex!(io::FilamentIO.VTKHDFFile, vs::AbstractVector{<:AbstractVector}, name::AbstractString)
Base.setindex!(io::FilamentIO.VTKHDFFile, data, name::AbstractString)
```

### Reading data

```@docs
read_vtkhdf
```

### Time series files

```@docs
TimeSeriesFile
Base.setindex!(tsf::TimeSeriesFile, filename::AbstractString, time::Real)
Base.empty!(tsf::TimeSeriesFile)
save(filename::AbstractString, tsf::TimeSeriesFile)
```

## Text format

The following functions can be used to dump filament positions (and other
quantities defined on filaments) to simple text files, and read positions back.
This can be a convenient way of making simulation checkpoints.

```@docs
write_to_text
read_from_text
```
