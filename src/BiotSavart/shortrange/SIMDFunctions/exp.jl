# Based on SLEEF.jl, adapted to allow using SIMD.Vec types.
# For scalars (or automatic vectorisation), performance is the same as directly using SLEEF.exp.

const MLN2E = 1.442695040888963407359924681001892137426645954152985934135449406931 # log2(e)

# Split log(2) into upper and lower parts
L2U(::Type{Float64}) = 0.69314718055966295651160180568695068359375
L2L(::Type{Float64}) = 0.28235290563031577122588448175013436025525412068e-12

L2U(::Type{Float32}) = 0.693145751953125f0
L2L(::Type{Float32}) = 1.428606765330187045f-06

max_exp(::Type{Float64}) = 709.78271114955742909217217426  # log 2^1023*(2-2^-52)
max_exp(::Type{Float32}) = 88.72283905206835f0             # log 2^127 *(2-2^-23)

min_exp(::Type{Float64}) = -7.451332191019412076235e2 # log 2^-1075
min_exp(::Type{Float32}) = -103.97208f0               # ≈ log 2^-150

@inline function exp_kernel(::Type{Float64}, x)
    c11 = 2.08860621107283687536341e-09
    c10 = 2.51112930892876518610661e-08
    c9  = 2.75573911234900471893338e-07
    c8  = 2.75572362911928827629423e-06
    c7  = 2.4801587159235472998791e-05
    c6  = 0.000198412698960509205564975
    c5  = 0.00138888888889774492207962
    c4  = 0.00833333333331652721664984
    c3  = 0.0416666666666665047591422
    c2  = 0.166666666666666851703837
    c1  = 0.50
    v = c11
    v = muladd(v, x, c10)
    v = muladd(v, x, c9)
    v = muladd(v, x, c8)
    v = muladd(v, x, c7)
    v = muladd(v, x, c6)
    v = muladd(v, x, c5)
    v = muladd(v, x, c4)
    v = muladd(v, x, c3)
    v = muladd(v, x, c2)
    v = muladd(v, x, c1)
    v
end

@inline function exp_kernel(::Type{Float32}, x)
    c6 = 0.000198527617612853646278381f0
    c5 = 0.00139304355252534151077271f0
    c4 = 0.00833336077630519866943359f0
    c3 = 0.0416664853692054748535156f0
    c2 = 0.166666671633720397949219f0
    c1 = 0.5f0
    v = c6
    v = muladd(v, x, c5)
    v = muladd(v, x, c4)
    v = muladd(v, x, c3)
    v = muladd(v, x, c2)
    v = muladd(v, x, c1)
    v
end

# Generalise unsafe_trunc to SIMD.Vec.
@generated function llvm_unsafe_trunc(::Type{I}, x::Vec{N, T}) where {I <: Integer, N, T <: AbstractFloat}
    I_llvm = SIMD.Intrinsics.llvm_type(I)
    T_llvm = SIMD.Intrinsics.llvm_type(T)
    s = """
    %2 = fptosi <$N x $T_llvm> %0 to <$N x $I_llvm>
    %3 = freeze <$N x $I_llvm> %2
    ret <$N x $I_llvm> %3
    """
    quote
        $(Expr(:meta, :inline))
        Vec{$N, $I}(Base.llvmcall($s, LVec{$N, $I}, Tuple{LVec{$N, $T}}, x.data))
    end
end

@inline _trunc(::Type{S}, x::AbstractFloat) where {S} = unsafe_trunc(S, x)

@inline function _trunc(::Type{S}, x::Vec{W, T}) where {S, W, T}
    llvm_unsafe_trunc(S, x)
end

@inline function ldexp2k(x, e::Union{Int, Vec{W, Int}}) where {W}
    T = ftype(x)
    a = pow2i(T, e >> 1)
    b = pow2i(T, e - (e >> 1))
    x = x * a
    x * b
end

@inline pow2i(::Type{T}, q::Int) where {T <: AbstractFloat} = integer2float(T, q + exponent_bias(T))
# @inline pow2i(::Type{T}, q::Vec{W, Int}) where {T <: AbstractFloat, W, Int} = convert(Vec{W, T}, one(UInt) << q)
@inline pow2i(::Type{T}, q::Vec) where {T <: AbstractFloat} = integer2float(T, q + exponent_bias(T))

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline function integer2float(::Type{Float64}, m::Vec{W, Int}) where {W}
    m_i64 = convert(Vec{W, Int64}, m)
    reinterpret(Vec{W, Float64}, m_i64 << significand_bits(Float64))
end

@inline function integer2float(::Type{Float32}, m::Vec{W, Int}) where {W}
    m_i32 = convert(Vec{W, Int32}, m)
    reinterpret(Vec{W, Float32}, m_i32 << significand_bits(Float32))
end

@inline exp(x) = _vexp(ftype(x), x)

@inline function _vexp(::Type{T}, d) where {T <: AbstractFloat}
    q = round(T(MLN2E) * d)
    qi = _trunc(Int, q)

    s = muladd(q, -L2U(T), d)
    s = muladd(q, -L2L(T), s)

    u = @inline exp_kernel(T, s)

    u = s * s * u + s + 1
    u = ldexp2k(u, qi)

    u = vifelse(d > max_exp(T), T(Inf), u)
    u = vifelse(d < min_exp(T), T(0), u)

    u
end
