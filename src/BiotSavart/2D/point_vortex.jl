export Vec2, PointVortices, ParamsPointVortex

const Vec2{T} = SVector{2, T}
const PointVortices{T} = AbstractVector{Vec2{T}}

include("params.jl")
