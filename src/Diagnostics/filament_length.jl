export filament_length

"""
    filament_length(f; quad = nothing) -> Real

Estimate length of one or more filaments.

Here `f` can be a [`Filament`](@ref) or a vector of filaments.

A quadrature rule may be optionally passed using `quad` (e.g. `quad = GaussLegendre(4)`),
which may result in better accuracy.
"""
function filament_length end

function filament_length(fs::AbstractVector{<:AbstractFilament}; kws...)
    T = eltype(eltype(eltype(fs)))
    @assert T <: AbstractFloat
    L = zero(T)
    for f ∈ fs
        L += filament_length(f; kws...)
    end
    L
end

function filament_length(f::AbstractFilament; quad = nothing)
    sum(segments(f)) do s
        Filaments.segment_length(s; quad)
    end
end
