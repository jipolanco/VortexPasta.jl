"""
    BiotSavart.SignValue

Enumerator (`@enum`) defining vortex signs `UP` and `DOWN`.

These are exported by the BiotSavart module, so one can directly use `UP` and `DOWN`.
"""
@enum SignValue::UInt8 begin
    UP   = 0x00
    DOWN = 0x01
end

@inline to_real(s::SignValue) = (s === UP) ? +1 : -1
@inline Base.:*(s::SignValue, x::Union{Number, SVector}) = to_real(s) * x  # automatically convert to ±1 when multiplying a number

function signs_to_charges!(charges::AbstractVector{T}, signs::AbstractVector{<:SignValue}) where {T <: Number}
    @assert eachindex(charges) == eachindex(signs)
    map!(charges, signs) do s
        @inline
        T(to_real(s))  # write result to `charges` (should work on GPU)
    end
end
