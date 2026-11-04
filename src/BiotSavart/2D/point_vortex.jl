export Vec2, PointVortices, ParamsPointVortex, PointVortexSign, SignValue, UP, DOWN

module PointVortexSign

"""
    BiotSavart.PointVortexSign.SignValue

Enumerator (`@enum`) defining vortex signs `UP` and `DOWN`.

These are exported by the BiotSavart module, so one can directly use `UP` and `DOWN`.
"""
@enum SignValue::UInt8 begin
    UP   = 0x00
    DOWN = 0x01
end

end  # module PointVortexSign

using .PointVortexSign: UP, DOWN, SignValue

const Vec2{T} = SVector{2, T}
const PointVortices{T} = AbstractVector{Vec2{T}}
const SignVector = AbstractVector{SignValue}

include("params.jl")
include("pointdata.jl")
include("cache.jl")
