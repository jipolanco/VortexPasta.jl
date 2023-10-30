# Implementation of cyclic pentadiagonal solver for periodic quintic spline interpolation.
# We use nearly pentadiagonal linear system algorithm by Lv and Le (Appl. Math and Comp.
# 2008), based on the Sherman–Morrison–Woodbury formula, which allows the use of standard
# pentadiagonal solvers (here we use a generic banded matrix solver for that).

using StaticArrays: SVector, SMatrix
using Bumper: Bumper, @no_escape, @alloc

# Specialisation for quintic splines.
function solve_spline_coefficients!(
        ::Val{6}, cs::PaddedVector{M}, ts::PaddedVector{M}, Xs::PaddedVector;
        Xoffset = zero(eltype(Xs)),
    ) where {M}
    periodise_coordinates!(cs, Xs, ts, Xoffset)  # useful if Xoffset ≠ 0
    _solve_quintic_spline_coefficients!(cs, ts)
    pad_periodic!(cs)
    cs
end

# Compute Uᵀ * y, where Uᵀ is a 2×m matrix which is only non-zero in its first and last
# 2 columns.
function lmul_utranspose(Uᵀ_left::SMatrix{2, 2}, Uᵀ_right::SMatrix{2, 2}, y::AbstractVector{<:Real})
    y_hi = SVector(y[begin], y[begin + 1])
    y_lo = SVector(y[end - 1], y[end])
    Uᵀ_left * y_hi + Uᵀ_right * y_lo
end

# Similar when `y` is a vector of SVector.
function lmul_utranspose(Uᵀ_left::SMatrix{2, 2}, Uᵀ_right::SMatrix{2, 2}, y::AbstractVector{<:SVector{N}}) where {N}
    y_hi = hcat(y[begin], y[begin + 1])' :: SMatrix{2, N}
    y_lo = hcat(y[end - 1], y[end])'     :: SMatrix{2, N}
    Uᵀ_left * y_hi + Uᵀ_right * y_lo
end

function _solve_quintic_spline_coefficients!(cs::AbstractVector{<:Vec3}, ts)
    n = length(ts)
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

function _solve_quintic_spline_coefficients!(buffers::NamedTuple, cs, ts)
    (; AB, V,) = buffers
    n = length(ts)
    m = n - 2
    kl, ku = 2, 2  # lower and upper band sizes

    fs = @view cs[begin:end - 2]
    fs_tilde = hcat(cs[end - 1], cs[end])' :: SMatrix{2, 3}

    (; A₂, Uᵀ_left, Uᵀ_right,) = _construct_quintic_spline_matrices!(AB, V, ts)

    Aband = SplineBandedMatrix(AB, m, kl, ku)
    banded_lu!(Aband)

    # Compute `r = f̃ - V A₂⁻¹ f`̃ and then solve `A y = r` in-place.
    gs = transpose(A₂ \ fs_tilde)
    @inbounds for i ∈ eachindex(fs)
        fs[i] = fs[i] - (V[i][1] * gs[:, 1] + V[i][2] * gs[:, 2])
    end

    # Now solve `A y = r` in-place. Also solve `A Z = V`.
    banded_ldiv!(Aband, fs, V)

    # Compute qⱼ = Uᵀzⱼ = Uᵀ A₁⁻¹ vⱼ
    Q = lmul_utranspose(Uᵀ_left, Uᵀ_right, V)
    B₂ = A₂ - Q  # = A₂ - Uᵀ A₁⁻¹ V

    # Compute x̂ solution onto y. This corresponds to coefficients 1:(n - 2).
    let y = fs
        local u = lmul_utranspose(Uᵀ_left, Uᵀ_right, y) :: SMatrix{2, 3}  # Uᵀy
        local w = B₂ \ u
        for i ∈ eachindex(y)
            @inbounds y[i] = y[i] + V[i][1] * w[1, :] + V[i][2] * w[2, :]
        end
    end

    # Compute x̃ solution. This corresponds to coefficients (n - 1):n.
    xs_tilde = let f = fs_tilde, y = fs
        local u = lmul_utranspose(Uᵀ_left, Uᵀ_right, y)  # Uᵀx̂
        local w = f - u
        A₂ \ w
    end

    # Finally, copy solution onto cs.
    cs[end - 1] = xs_tilde[1, :]
    cs[end - 0] = xs_tilde[2, :]

    cs
end

@inbounds function _construct_quintic_spline_matrices!(
        A::AbstractMatrix{T}, V::AbstractVector{<:SVector{2, T}}, ts,
    ) where {T}
    @assert size(A, 1) == 5  # = kl + ku + 1
    D = 3  # = ku + 1

    m = size(A, 2)
    n = m + 2

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
            A[D + i - j, j] = bs[4 - j]
        end
    end

    # Second row
    let i = 2
        bs = eval_bsplines(order, ts, i)
        V[i] = (zero(T), bs[5])
        for j ∈ 1:4
            A[D + i - j, j] = bs[5 - j]
        end
    end

    # Rows 3:(n - 4)
    for i ∈ 3:(n - 4)
        bs = eval_bsplines(order, ts, i)
        for l ∈ eachindex(bs)
            j = i + 3 - l
            A[D - 3 + l, j] = bs[l]
        end
    end

    # Row n - 3
    let i = n - 3
        bs = eval_bsplines(order, ts, i)
        V[i] = (bs[1], zero(T))
        for l ∈ 2:5
            j = i + 3 - l
            A[D - 3 + l, j] = bs[l]
        end
    end

    # Row n - 2
    let i = n - 2
        bs = eval_bsplines(order, ts, i)
        V[i] = (bs[2], bs[1])
        for l ∈ 3:5
            j = i + 3 - l
            A[D - 3 + l, j] = bs[l]
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
