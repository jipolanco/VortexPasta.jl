# Implementation of cyclic pentadiagonal solver for periodic quintic spline interpolation.
# We use nearly pentadiagonal linear system algorithm by Lv and Le (Appl. Math and Comp.
# 2008), based on the Sherman–Morrison–Woodbury formula, which allows the use of standard
# pentadiagonal solvers (here we use a generic banded matrix solver for that).

using StaticArrays: SVector, SMatrix
using Bumper: Bumper, @no_escape, @alloc

# Specialisation for quintic splines.
function solve_spline_coefficients!(
        ::Val{6}, cs::PaddedVector{M}, ts::PaddedVector{M}, Xs::AbstractVector;
        Xoffset = zero(eltype(Xs)),
    ) where {M}
    periodise_coordinates!(cs, Xs, ts, Xoffset)  # useful if Xoffset ≠ 0
    _solve_quintic_spline_coefficients!(cs, ts)
    pad_periodic!(cs)
    cs
end

# Compute Uᵀ * y, where Uᵀ is a 2×m matrix which is only non-zero in its first and last
# 2 columns. Here `y` is a vector of SVector.
function lmul_utranspose(Uᵀ_left::SMatrix{2, 2}, Uᵀ_right::SMatrix{2, 2}, y::AbstractVector)
    y_hi = hcat_transpose(y[begin], y[begin + 1])
    y_lo = hcat_transpose(y[end - 1], y[end])
    Uᵀ_left * y_hi + Uᵀ_right * y_lo
end

function _solve_quintic_spline_coefficients!(cs::AbstractVector, ts)
    n = length(ts)
    if n < 7
        _solve_quintic_spline_coefficients_small!(cs, ts)
        return cs
    end
    m = n - 2
    kl, ku = 2, 2  # lower and upper band sizes
    T = eltype(ts)
    buf = Bumper.default_buffer()
    @no_escape buf begin
        V = @alloc(SVector{2, T}, m)
        AB = @alloc(T, kl + ku + 1, m)
        buffers = (; AB, V,)
        _solve_quintic_spline_coefficients!(buffers, cs, ts)
    end
    cs
end

hcat_transpose(u::SVector{N}, v::SVector{N}) where {N} = hcat(u, v)' :: SMatrix{2, N}
hcat_transpose(u::T, v::T) where {T <: Number} = SVector{2}(u, v)  # scalar case

quintic_ldiv_terms(ys::SMatrix{2}) = (ys[1, :], ys[2, :])
quintic_ldiv_terms(ys::SVector{2}) = ys

function _solve_quintic_spline_coefficients!(buffers::NamedTuple, cs::AbstractVector, ts)
    (; AB, V,) = buffers
    n = length(ts)
    m = n - 2
    kl, ku = 2, 2  # lower and upper band sizes

    fs = @view cs[begin:end - 2]
    fs_tilde = hcat_transpose(cs[end - 1], cs[end])

    (; A₂, Uᵀ_left, Uᵀ_right,) = _construct_quintic_spline_matrices!(AB, V, ts)

    Aband = SplineBandedMatrix(AB, m, kl, ku)
    banded_lu!(Aband)

    # Compute `r = f̃ - V A₂⁻¹ f`̃ and then solve `A y = r` in-place.
    gs = quintic_ldiv_terms(A₂ \ fs_tilde)
    @inbounds for i ∈ eachindex(fs)
        fs[i] = fs[i] - (V[i][1] * gs[1] + V[i][2] * gs[2])
    end

    # Now solve `A y = r` in-place. Also solve `A Z = V`.
    banded_ldiv!(Aband, fs, V)

    # Compute qⱼ = Uᵀzⱼ = Uᵀ A₁⁻¹ vⱼ
    Q = lmul_utranspose(Uᵀ_left, Uᵀ_right, V)
    B₂ = A₂ - Q  # = A₂ - Uᵀ A₁⁻¹ V

    # Compute x̂ solution onto y. This corresponds to coefficients 1:(n - 2).
    let y = fs
        local u = lmul_utranspose(Uᵀ_left, Uᵀ_right, y)  # Uᵀy
        local ws = quintic_ldiv_terms(B₂ \ u)
        for i ∈ eachindex(y)
            @inbounds y[i] = y[i] + V[i][1] * ws[1] + V[i][2] * ws[2]
        end
    end

    # Compute x̃ solution. This corresponds to coefficients (n - 1):n.
    xs_tilde = let f = fs_tilde, y = fs
        local u = lmul_utranspose(Uᵀ_left, Uᵀ_right, y)  # Uᵀx̂
        local w = f - u
        quintic_ldiv_terms(A₂ \ w)
    end

    # Finally, copy solution onto cs.
    cs[end - 1] = xs_tilde[1]
    cs[end - 0] = xs_tilde[2]

    cs
end

function _construct_quintic_spline_matrices!(
        A::AbstractMatrix{T}, V::AbstractVector{<:SVector{2, T}}, ts,
    ) where {T}
    @assert size(A, 1) == 5  # = kl + ku + 1
    D = 3  # = ku + 1

    m = size(A, 2)
    @assert m ≥ 5  # this is ensured by minimum_nodes(::QuinticSplineMethod) = 7
    n = m + 2      # this is the number of nodes

    # Set arrays to zero
    fill!(A, zero(T))
    fill!(V, zero(eltype(V)))

    order = Val(6)

    # First row
    let i = 1
        # Note: B-splines are returned in reverse order!!
        # In this case: (b[5], b[4], ..., b[1]).
        bs = eval_bsplines(order, ts, i) :: NTuple{5}
        V[i] = (bs[5], bs[4])
        for j ∈ 1:3
            @inbounds A[D + i - j, j] = bs[4 - j]
        end
    end

    # Second row
    let i = 2
        bs = eval_bsplines(order, ts, i)
        V[i] = (zero(T), bs[5])
        for j ∈ 1:4
            @inbounds A[D + i - j, j] = bs[5 - j]
        end
    end

    # Rows 3:(n - 4)
    for i ∈ 3:(n - 4)
        bs = eval_bsplines(order, ts, i)
        for l ∈ eachindex(bs)
            j = i + 3 - l
            @inbounds A[D - 3 + l, j] = bs[l]
        end
    end

    # Row n - 3
    let i = n - 3
        bs = eval_bsplines(order, ts, i)
        V[i] = (bs[1], zero(T))
        for l ∈ 2:5
            j = i + 3 - l
            @inbounds A[D - 3 + l, j] = bs[l]
        end
    end

    # Row n - 2
    let i = n - 2
        bs = eval_bsplines(order, ts, i)
        V[i] = (bs[2], bs[1])
        for l ∈ 3:5
            j = i + 3 - l
            @inbounds A[D - 3 + l, j] = bs[l]
        end
    end

    # Rows (n - 1):n
    bsₙ₋₁ = eval_bsplines(order, ts, n - 1)
    bsₙ = eval_bsplines(order, ts, n)

    A₂ = SMatrix{2, 2}(
        bsₙ₋₁[3], bsₙ[4],  # first column
        bsₙ₋₁[2], bsₙ[3],  # second column
    )

    Uᵀ_left = SMatrix{2, 2}(
        bsₙ₋₁[1], bsₙ[2],  # first column
               0, bsₙ[1],  # second column
    )

    Uᵀ_right = SMatrix{2, 2}(
        bsₙ₋₁[5],      0,  # first column
        bsₙ₋₁[4], bsₙ[5],  # second column
    )

    (; A₂, Uᵀ_left, Uᵀ_right,)
end

## ================================================================================ ##

# Special case of small linear systems (n = 5, 6), which are not supported by the above
# solver.
# We directly solve the linear system. Note that a cyclic pentadiagonal matrix is basically
# a full matrix for these sizes, so we can use generic solvers.

function _solve_quintic_spline_coefficients_small!(cs, ts)
    n = length(cs)
    @assert 5 ≤ n ≤ 6
    if n == 6
        _solve_quintic_spline_coefficients_small!(cs, ts, Val(6))
    elseif n == 5
        _solve_quintic_spline_coefficients_small!(cs, ts, Val(5))
    end
    nothing
end

function _solve_quintic_spline_coefficients_small!(cs, ts, ::Val{n}) where {n}
    A = _construct_quintic_spline_matrix_small(ts, Val(n))
    ys = SVector{n}(cs)
    xs = A \ ys
    copyto!(cs, xs)
end

function _construct_quintic_spline_matrix_small(ts, ::Val{n}) where {n}
    @assert n == length(ts)
    @assert n ≥ 5
    n² = n * n
    T = eltype(ts)
    A = zero(MMatrix{n, n, T, n²})
    order = Val(6)
    for i ∈ 1:n
        bs = eval_bsplines(order, ts, i) :: NTuple{5}
        for (ib, b) ∈ pairs(bs)
            j = i + 3 - ib
            if j < 1
                j += n
            elseif j > n
                j -= n
            end
            A[i, j] = b
        end
    end
    SMatrix(A)
end
