# Banded matrix solvers for spline interpolation matrices.
# Based on Carl de Boor's routines (BANSLV, BANFAC) originally in FORTRAN77.

using LinearAlgebra: ZeroPivotException

struct SplineBandedMatrix{A <: AbstractArray}
    w      :: A    # data
    nrow   :: Int
    nbandl :: Int
    nbandu :: Int
end

# Perform in-place LU factorisation of collocation matrix without pivoting.
#
# Takes advantage of the totally positive property of collocation matrices
# appearing in spline calculations (de Boor 1978).
#
# The code is ported from Carl de Boor's BANFAC routine in FORTRAN77, via its
# [FORTRAN90 version by John Burkardt](https://people.math.sc.edu/Burkardt/f_src/pppack/pppack.html).
#
# Adapted from version in BSplineKit.jl.
# Simplified by removing triangular matrix cases (either with nbandl = 0 or nbandu = 0).
function banded_lu!(A::SplineBandedMatrix)
    (; w, nrow, nbandl, nbandu,) = A
    nroww = size(w, 1)
    @assert nrow == size(w, 2) "expected a square banded matrix"
    @assert nroww == nbandl + nbandu + 1 "inconsistent number of bands"
    @assert nbandl > 0 && nbandu > 0 "expected both lower and upper bands"
    middle = nbandu + 1  # w[middle, :] contains the main diagonal of A

    if nrow == 1
        @inbounds iszero(w[middle, nrow]) && throw(ZeroPivotException(1))
        return A
    end

    @inbounds for i = 1:(nrow - 1)
        pivot = w[middle, i]  # pivot for the i-th step
        iszero(pivot) && throw(ZeroPivotException(i))
        ipiv = inv(pivot)

        # Divide each entry in column `i` below the diagonal by `pivot`.
        for j = 1:min(nbandl, nrow - i)
            w[middle + j, i] *= ipiv
        end

        # Subtract A[i, i+k] * (i-th column) from (i+k)-th column (below row `i`).
        for k = 1:min(nbandu, nrow - i)
            factor = w[middle - k, i + k]
            for j = 1:min(nbandl, nrow - i)
                w[middle - k + j, i + k] -= factor * w[middle + j, i]
            end
        end
    end

    # Check the last diagonal entry.
    @inbounds iszero(w[middle, nrow]) && throw(ZeroPivotException(nrow))

    A
end

# Solution of banded linear system A * x = y.
# The code is adapted from Carl de Boor's BANSLV routine in FORTRAN77, via its
# FORTRAN90 version by John Burkardt.
#
# Here `A` actually contains the LU decomposition of the original banded matrix (obtained
# from `banded_lu!`).
#
# We also allow multiple RHSs. We use a @generated function to ensure type stability
# when the RHSs have different element types (which is actually our case).
@generated function banded_ldiv!(A::SplineBandedMatrix, xs::Vararg{AbstractVector})
    ex_z1 = :()
    ex_zi = :()
    ex_fw = :()
    ex_bw = :()
    for k ∈ eachindex(xs)
        ex_z1 = quote
            $ex_z1
            xs[$k][1] /= z
        end
        ex_zi = quote
            $ex_zi
            xs[$k][i] /= z
        end
        ex_fw = quote
            $ex_fw
            xs[$k][i + j] -= xs[$k][i] * y
        end
        ex_bw = quote
            $ex_bw
            xs[$k][i - j] -= xs[$k][i] * y
        end
    end

    ex = quote
        (; w, nrow, nbandl, nbandu,) = A
        middle = nbandu + 1

        if nrow == 1
            let z = w[middle, 1]
                $ex_z1
            end
            return xs
        end

        # Forward pass:
        #
        # For i = 1:(nrow-1), subtract RHS[i] * (i-th column of L) from the
        # right hand side, below the i-th row.
        if nbandl != 0
            @inbounds for i = 1:(nrow - 1)
                jmax = min(nbandl, nrow - i)
                for j = 1:jmax
                    let y = w[middle + j, i]
                        $ex_fw
                    end
                end
            end
        end

        # Backward pass:
        #
        # For i = nrow:-1:1, divide RHS[i] by the i-th diagonal entry of
        # U, then subtract RHS[i]*(i-th column of U) from right hand side, above the
        # i-th row.
        @inbounds for i = nrow:-1:2
            let z = w[middle, i]
                $ex_zi
            end
            for j = 1:min(nbandu, i - 1)
                let y = w[middle - j, i]
                    $ex_bw
                end
            end
        end

        let z = w[middle, 1]
            $ex_z1
        end

        xs
    end
end
