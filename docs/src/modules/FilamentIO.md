# FilamentIO

```@meta
CurrentModule = VortexPasta.FilamentIO
```

```@docs
FilamentIO
```

## Writing data

```@docs
write_vtkhdf
Base.setindex!(io::VTKHDFFile, vs::AbstractVector{<:AbstractVector}, name::AbstractString)
Base.setindex!(io::VTKHDFFile, data, name::AbstractString)
```

## Reading data

```@docs
read_vtkhdf
```

## Index

```@index
Pages = ["FilamentIO.md"]
```
