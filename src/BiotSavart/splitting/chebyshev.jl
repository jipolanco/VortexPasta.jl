# Approximate smooth functions using Chebyshev series.

module ChebyshevApproximations

using Adapt: Adapt
using FFTW: FFTW

# Approximates a function in interval [-L, L] using Chebyshev series.
# It can take symmetries into account (odd/even) to reduce the number of operations.
# In that case, it is *assumed* that the zero coefficients have already been discarded
# before constructing the series.
struct ChebyshevSeries{Symmetry, T, Coefs <: AbstractVector{T}} <: Function
    cs::Coefs
    L::T
    Linv::T  # = 1/L
    function ChebyshevSeries{Symmetry}(cs, L, Linv) where {Symmetry}
        T = eltype(cs)
        new{Symmetry, T, typeof(cs)}(cs, T(L), T(Linv))
    end
end

function ChebyshevSeries(cs, L; symmetry::Val{Symmetry} = Val(:none)) where {Symmetry}
    ChebyshevSeries{Symmetry}(cs, L, 1 / L)
end

@inline function Adapt.adapt_structure(to, g::ChebyshevSeries{Symmetry}) where {Symmetry}
    (; L, Linv) = g
    cs = Adapt.adapt(to, g.cs)
    ChebyshevSeries{Symmetry}(cs, L, Linv)
end

@inline (f::ChebyshevSeries{:none})(x) = eval_chebyshev_clenshaw(f.cs, x * f.Linv)
@inline (f::ChebyshevSeries{:even})(x) = eval_chebyshev_clenshaw_even(f.cs, x * f.Linv)
@inline (f::ChebyshevSeries{:odd})(x) = eval_chebyshev_clenshaw_odd(f.cs, x * f.Linv)

@inline function eval_chebyshev_clenshaw(cs::AbstractVector, x)
    Base.require_one_based_indexing(cs)
    N = length(cs)
    N == 0 && return zero(x)
    Bpp = zero(x)
    Bp = zero(x)
    twox = 2 * x
    for n in N:-1:2
        @inbounds B = muladd(twox, Bp, cs[n] - Bpp)
        Bpp, Bp = Bp, B
    end
    @inbounds cs[begin] + x * Bp - Bpp
end

@inline function eval_chebyshev_clenshaw_even(cs::AbstractVector, x)
    Base.require_one_based_indexing(cs)
    N = length(cs)
    N == 0 && return zero(x)
    Bpp = zero(x)
    Bp = zero(x)
    twox = 2 * x
    @inbounds for n in N:-1:2
        B = muladd(twox, Bp, cs[n] - Bpp)
        Bpp, Bp = Bp, B
        B = muladd(twox, Bp, -Bpp)
        Bpp, Bp = Bp, B
    end
    @inbounds cs[begin] + x * Bp - Bpp
end

@inline function eval_chebyshev_clenshaw_odd(cs::AbstractVector, x)
    Base.require_one_based_indexing(cs)
    N = length(cs)
    N == 0 && return zero(x)
    Bpp = zero(x)
    Bp = @inbounds one(x) * cs[N]
    twox = 2 * x
    @inbounds for n in (N - 1):-1:1
        B = muladd(twox, Bp, -Bpp)
        Bpp, Bp = Bp, B
        B = muladd(twox, Bp, cs[n] - Bpp)
        Bpp, Bp = Bp, B
    end
    x * Bp - Bpp
end

function integrate(f::ChebyshevSeries{:none})
    (; cs, L) = f
    Base.require_one_based_indexing(cs)
    Nint = length(cs) + 1
    cs_int = similar(cs, Nint)
    δ = 1  # for zero -> one-based indexing
    cs_int[δ + 0] = 0
    cs_int[δ + 1] = (2 * cs[δ + 0] - cs[δ + 2]) / 2
    for n in 2:(Nint - 3)
        cs_int[δ + n] = (cs[δ + n - 1] - cs[δ + n + 1]) / (2 * n)
    end
    let n = Nint - 2
        cs_int[δ + n] = cs[δ + n - 1] / (2 * n)
    end
    let n = Nint - 1
        cs_int[δ + n] = cs[δ + n - 1] / (2 * n)
    end
    cs_int .*= L  # rescale according to interval size
    ChebyshevSeries(cs_int, L; symmetry = Val(:none))
end

function integrate(f::ChebyshevSeries{:even})
    (; cs, L) = f
    Base.require_one_based_indexing(cs)
    N = length(cs)
    cs_int = similar(cs, N)
    cs_int[1] = (2 * cs[1] - cs[2]) / 2  # T[1] <- T[0], T[2]
    for i in 2:(N - 1)
        n = 2i - 1
        # i = 2: T[3] <- T[2], T[4]
        cs_int[i] = (cs[i] - cs[i + 1]) / (2 * n)
    end
    let i = N
        n = 2i - 1
        cs_int[i] = cs[i] / (2 * n)
    end
    cs_int .*= L  # rescale according to interval size
    ChebyshevSeries(cs_int, L; symmetry = Val(:odd))
end

# Approximate function f(x) using Chebyshev series.
# Here N is the initial number of coefficients.
# It will be increased adaptatively if the requested precision is not achieved.
function approximate(
        f::F, L::T;
        symmetry::Val{Symmetry} = Val(:none), N = 64, rtol = 10 * eps(T),
    ) where {F <: Function, T <: AbstractFloat, Symmetry}
    # Using extrema grid
    xs = [L * cospi(n / (N - 1)) for n in 0:(N - 1)]
    us = f.(xs)
    cs = FFTW.r2r(us, FFTW.REDFT00)
    cs[begin] /= 2
    cs[end] /= 2
    cs ./= (N - 1)

    # Using roots grid
    # xs = [L * cospi((n - 1/2) / N) for n in 1:N]
    # cs = FFTW.r2r(us_eval, FFTW.REDFT10)
    # cs[begin] /= 2
    # cs ./= N

    # Determine where to truncate the series.
    # Note that, if the original function has even or odd symmetry, then half the
    # coefficients will be practically zero, so we need to check at least two consecutive
    # coefficients.
    cs_sum = sum(abs, cs)
    atol = rtol * cs_sum
    n = firstindex(cs)
    while n < lastindex(cs)
        if abs(cs[n]) < atol && abs(cs[n + 1]) < atol  # if two consecutive coefs are zero
            break
        end
        n += 1
    end

    if n == lastindex(cs)
        # We need more coefficients, try again with larger N.
        return approximate(f, L; N = 2 * N, rtol, symmetry)
    end

    @assert n > firstindex(cs)
    cs_final = if Symmetry === :even
        cs[1:2:(n - 1)]
    elseif Symmetry === :odd
        cs[2:2:(n - 1)]
    else
        cs[1:1:(n - 1)]
    end

    ChebyshevSeries(cs_final, L; symmetry)
end

end
