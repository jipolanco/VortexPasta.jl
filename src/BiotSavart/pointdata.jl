# This stores the locations s⃗ and charges q * s⃗′ used to compute both short-range and
# long-range interactions, allowing to reuse computations. Note that locations and charges
# are obtained via interpolation in-between filament nodes.
#
# This is also reused by long-range interactions to perform interpolations from Fourier to
# physical space (see `interpolate_to_physical!`).
# In that case, `points` contains the interpolation points (usually the filament nodes) and
# `charges` the interpolation values (usually velocities or streamfunction values).
struct PointData{
        T <: AbstractFloat,
        Points <: StructVector{Vec3{T}},
        # Note: complex may be needed by long-range backend (even though values are always real!)
        Charges <: StructVector{Vec3{Complex{T}}},
    }
    points  :: Points   # interpolated locations s⃗ on segments
    charges :: Charges  # rescaled tangent vector q * s⃗′ on segments (where `q` is the quadrature weight)
end

function PointData(::Type{T}) where {T}
    points = StructVector{Vec3{T}}(undef, 0)
    charges = StructVector{Vec3{Complex{T}}}(undef, 0)
    PointData(points, charges)
end

"""
    set_num_points!(data::PointData, N::Integer)

Set the total number of non-uniform points that the cache must hold.

This will reallocate space to make all points fit in the cache. It will also reset the
contributions of previously-added charges.

Must be called before [`add_point!`](@ref).
"""
function set_num_points!(data::PointData, N)
    resize!(data.points, N)
    resize!(data.charges, N)
    data
end

function _count_charges(quad::AbstractQuadrature, fs::AbstractVector{<:ClosedFilament})
    Nq = length(quad)        # number of evaluation points per filament segment
    Np = sum(f -> length(segments(f)), fs)  # total number of segments among all filaments (assumes closed filaments!!)
    Np * Nq
end

"""
    add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, quad::AbstractQuadrature)

Add vector charges at multiple non-uniform locations.

This is used in particular for type-1 NUFFTs, to transform from non-uniform data in physical
space to uniform data in Fourier space. It must be called before [`compute_vorticity_fourier!`](@ref).
"""
function add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, quad::AbstractQuadrature)
    Ncharges = _count_charges(quad, fs)
    set_num_points!(data, Ncharges)
    n = 0
    for f ∈ fs
        n = _add_point_charges!(data, f, n, quad)
    end
    @assert n == Ncharges
    nothing
end

function _add_point_charges!(data::PointData, f, n, quad::AbstractQuadrature)
    ζs, ws = quadrature(quad)
    ts = knots(f)
    @inbounds for i ∈ eachindex(segments(f))
        Δt = ts[i + 1] - ts[i]
        for (ζ, w) ∈ zip(ζs, ws)
            s⃗ = f(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
            # Note: the vortex circulation Γ is included in the Ewald operator and
            # doesn't need to be included here.
            q = w * Δt
            add_pointcharge!(data, s⃗, q * s⃗′, n += 1)
        end
    end
    n
end

function _add_point_charges!(data::PointData, f, n, ::NoQuadrature)
    @inbounds for i ∈ eachindex(segments(f))
        s⃗ = (f[i] + f[i + 1]) ./ 2
        s⃗′_dt = f[i + 1] - f[i]
        add_pointcharge!(data, s⃗, s⃗′_dt, n += 1)
    end
    n
end

function add_pointcharge!(data::PointData, X::Vec3, Q::Vec3, i::Int)
    @inbounds data.points[i] = X
    @inbounds data.charges[i] = Q
    data
end

"""
    add_point!(data::PointData, X::Vec3, i::Int)

Add an interpolation point for type-2 NUFFT.

This is used for type-2 NUFFTs, to transform (interpolate) from uniform data in Fourier
space to non-uniform data in physical space.

The total number of locations must be first set via [`set_num_points!`](@ref).
"""
function add_point!(data::PointData, X::Vec3, i::Int)
    @inbounds data.points[i] = X
    data
end
