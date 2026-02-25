# Note: all input values may be SIMD types (Vec or tuple of Vec).
@inline function short_range_integrand(::Velocity, a, b, r_inv::V, qs⃗′::NTuple{3, V}, r⃗::NTuple{3, V}) where {V}
    factor = (a + b) * r_inv^3
    vec = crossprod(qs⃗′, r⃗)
    map(vec) do component
        @inline
        factor * component
    end
end

@inline function short_range_integrand(::Streamfunction, a, b, r_inv, qs⃗′, r⃗)
    factor = a * r_inv
    vec = qs⃗′
    map(vec) do component
        @inline
        factor * component
    end
end

@inline function long_range_integrand(::Velocity, a, b, r_inv::V, qs⃗′::NTuple{3, V}, r⃗::NTuple{3, V}) where {V}
    factor = (a - b) * r_inv^3
    vec = crossprod(qs⃗′, r⃗)
    map(vec) do component
        @inline
        factor * component
    end
end

@inline function long_range_integrand(::Streamfunction, a, b, r_inv, qs⃗′, r⃗)
    factor = a * r_inv
    vec = qs⃗′
    map(vec) do component
        @inline
        factor * component
    end
end

@inline function full_integrand(::Velocity, r_inv, qs⃗′::NTuple{3, V}, r⃗::NTuple{3, V}) where {V}
    factor = r_inv^3
    vec = crossprod(qs⃗′, r⃗)
    map(vec) do component
        @inline
        factor * component
    end
end

@inline function full_integrand(::Streamfunction, r_inv, qs⃗′::NTuple{3, V}, r⃗::NTuple{3, V}) where {V}
    factor = r_inv
    vec = qs⃗′
    map(vec) do component
        @inline
        factor * component
    end
end

# Note: V may be a real or a SIMD vector.
@inline function crossprod(u::T, v::T) where {T <: NTuple{3}}
    @inbounds (
        u[2] * v[3] - u[3] * v[2],
        u[3] * v[1] - u[1] * v[3],
        u[1] * v[2] - u[2] * v[1],
    ) :: T
end
