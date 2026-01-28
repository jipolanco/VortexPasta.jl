module SIMDFunctions

using SIMD: SIMD, Vec, LVec, vifelse
using Base.Math: exponent_bias, significand_bits

public erf, exp

function erf end
function exp end

include("erf.jl")
include("exp.jl")

end
