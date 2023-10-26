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

function _solve_quintic_spline_coefficients!(cs::AbstractVector{<:Vec3}, ts)
    n = length(ts)
    m = n - 2
    T = eltype(ts)
    buf = Bumper.default_buffer()

    # Compute Uᵀ * y, where Uᵀ is a 2×m matrix which is only non-zero in its first and last
    # 2 columns.
    function lmul_utranspose(Uᵀ_left::SMatrix{2, 2}, Uᵀ_right::SMatrix{2, 2}, y)
        y_hi = SVector(y[begin], y[begin + 1])
        y_lo = SVector(y[end - 1], y[end])
        Uᵀ_left * y_hi + Uᵀ_right * y_lo
    end

    @no_escape buf begin
        vs = ntuple(i -> @alloc(T, m), Val(2))
        fs = ntuple(i -> @alloc(T, m), Val(3))

        kl, ku = 2, 2  # lower and upper band sizes
        AB = @alloc(T, kl + ku + 1, m)
        (; A₂, Uᵀ_left, Uᵀ_right,) = _construct_quintic_spline_matrices!(AB, vs, ts)
        fs_tilde = _construct_quintic_spline_rhs!(fs, cs)

        Aband = SplineBandedMatrix(AB, m, kl, ku)
        banded_lu!(Aband)

        # Compute `r = f̃ - V A₂⁻¹ f`̃ and then solve `A y = r` in-place.
        for (f, f_tilde) ∈ zip(fs, fs_tilde)
            g = A₂ \ f_tilde  # = A₂⁻¹ f̃
            @inbounds for i ∈ eachindex(f)
                f[i] = f[i] - (vs[1][i] * g[1] + vs[2][i] * g[2])
            end
            banded_ldiv!(Aband, f)
        end
        ys = fs  # define alias for convenience

        # Also solve `A zⱼ = vⱼ` in-place.
        for v ∈ vs
            banded_ldiv!(Aband, v)
        end
        zs = vs  # define alias for convenience

        # Compute qⱼ = Uᵀzⱼ = Uᵀ A₁⁻¹ vⱼ
        qs = map(zs) do y
            lmul_utranspose(Uᵀ_left, Uᵀ_right, y)
        end
        B₂ = A₂ - hcat(qs...)::SMatrix{2, 2}  # = A₂ - Uᵀ A₁⁻¹ V

        # Compute x̂ solution onto y. This corresponds to coefficients 1:(n - 2).
        for y ∈ ys
            local u = lmul_utranspose(Uᵀ_left, Uᵀ_right, y)  # Uᵀy
            local w = B₂ \ u
            for i ∈ eachindex(y)
                @inbounds y[i] = y[i] + zs[1][i] * w[1] + zs[2][i] * w[2]
            end
        end

        # Compute x̃ solution. This corresponds to coefficients (n - 1):n.
        xs_tilde = map(fs_tilde, ys) do f, y
            local u = lmul_utranspose(Uᵀ_left, Uᵀ_right, y)  # Uᵀx̂
            local w = f - u
            A₂ \ w
        end

        # Finally, copy solution onto cs.
        _copy_quintic_spline_solution!(cs, ys, xs_tilde)
    end

    cs
end

function _copy_quintic_spline_solution!(cs::AbstractVector{<:Vec3}, xs::NTuple{3}, xs_tilde::NTuple{3})
    @assert length(xs[1]) == length(cs) - 2
    for i ∈ eachindex(cs)[begin:end - 2]
        cs[i] = map(x -> x[i], xs)
    end
    cs[end - 1] = map(x -> x[1], xs_tilde)
    cs[end - 0] = map(x -> x[2], xs_tilde)
    cs
end

@inbounds function _construct_quintic_spline_rhs!(fs::NTuple{3}, cs::AbstractVector{<:Vec3})
    @assert length(fs[1]) == length(cs) - 2
    for i ∈ eachindex(cs)[begin:end - 2]
        for (j, cj) ∈ pairs(cs[i])
            fs[j][i] = cj
        end
    end
    fs_tilde = ntuple(Val(3)) do j
        SVector(cs[end - 1][j], cs[end][j])
    end
    fs_tilde
end

@inbounds function _construct_quintic_spline_matrices!(
        A::AbstractMatrix, vs::NTuple{2}, ts,
    )
    @assert size(A, 1) == 5  # = kl + ku + 1
    D = 3  # = ku + 1

    T = eltype(A)
    m = size(A, 2)
    n = m + 2

    # Set arrays to zero
    fill!(A, zero(T))
    for v ∈ vs
        fill!(v, zero(T))
    end

    order = Val(6)

    # First row
    let i = 1
        # Note: B-splines are returned in reverse order!!
        # In this case: (b[5], b[4], ..., b[1]).
        bs = eval_bsplines(order, ts, i) :: NTuple{5}
        vs[1][i] = bs[5]
        vs[2][i] = bs[4]
        for j ∈ 1:3
            A[D + i - j, j] = bs[4 - j]
        end
    end

    # Second row
    let i = 2
        bs = eval_bsplines(order, ts, i)
        vs[2][i] = bs[5]
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
        vs[1][i] = bs[1]
        for l ∈ 2:5
            j = i + 3 - l
            A[D - 3 + l, j] = bs[l]
        end
    end

    # Row n - 2
    let i = n - 2
        bs = eval_bsplines(order, ts, i)
        vs[2][i] = bs[1]
        vs[1][i] = bs[2]
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
