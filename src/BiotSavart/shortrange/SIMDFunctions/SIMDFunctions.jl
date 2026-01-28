module SIMDFunctions

using SIMD: SIMD, Vec, LVec, vifelse
using Base.Math: exponent_bias, significand_bits

public exp, erf, erfc

function exp end
function erf end
function erfc end

ftype(::T) where {T <: AbstractFloat} = T
ftype(::Vec{N, T}) where {N, T <: AbstractFloat} = T

include("exp.jl")
include("erf.jl")
include("erfc.jl")

end
