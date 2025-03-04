export write_to_text, read_from_text

"""
    write_to_text(io::IO, fs::AbstractVector{<:AbstractVector})
    write_to_text(filename::AbstractString, fs::AbstractVector{<:AbstractVector})

Dump filament locations or other filament data to text file.

Here `fs` can be a list of filaments or a list of values (e.g. velocities) on filament nodes.

The output is a text file where:

- each line represents a single filament,

- node coordinates (or vector values such as velocity) are flattened so that they can be
  written in a single line,

- the end node is included. For a closed filament, this means that the first and last written points are identical.
  This also means that we write `N + 1` values for a filament of length `N` (= number of _independent_ nodes).

One can use [`read_from_text`](@ref) to generate filaments from files written by this function.
"""
function write_to_text(io::IO, fs::AbstractVector{<:AbstractVector}; delim = '\t')
    for f in fs
        for i in 0:length(f)
            xs = f[begin + i]
            for x in xs
                print(io, x, delim)
            end
        end
        print(io, '\n')
    end
    nothing
end

write_to_text(filename::AbstractString, fs::AbstractVector{<:AbstractVector}; kws...) =
    open(io -> write_to_text(io, fs; kws...), filename, "w")

"""
    read_from_text(io::IO, ::Type{T}, method::DiscretisationMethod) -> fs
    read_from_text(filename::AbstractString, ::Type{T}, method::DiscretisationMethod) -> fs

Read filament locations from text file.

This function can be used to read filament locations written by [`write_to_text`](@ref).

The type `T` corresponds to the wanted numerical precision. It is usually `Float64` or `Float32`.
"""
function read_from_text(io::IO, ::Type{T}, method::DiscretisationMethod; delim = '\t') where {T}
    FilamentType = typeof(Filaments.init(ClosedFilament{T}, 0, method))
    fs = FilamentType[]
    data = Vector{T}(undef, 0)
    while !eof(io)
        s = readuntil(io, delim)
        x = parse(T, s)
        push!(data, x)
        if Char(peek(io)) == '\n'  # end of filament
            points = reinterpret(Vec3{T}, data)
            offset = points[end] - points[begin]
            f = Filaments.init(ClosedFilament, @view(points[begin:end-1]), method; offset)
            push!(fs, f)
            empty!(data)
            read(io, Char)  # == '\n'
        end
    end
    fs
end

read_from_text(filename::AbstractString, ::Type{T}, method::DiscretisationMethod; kws...) where {T} =
    open(io -> read_from_text(io, T, method; kws...), filename, "r")
