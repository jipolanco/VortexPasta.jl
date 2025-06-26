module VortexPastaMakieExt

using VortexPasta.Filaments
import VortexPasta.Filaments: filamentplot, filamentplot!  # make sure @recipe overloads functions in VP.Filaments

using Makie: Makie, Attributes, Observable, @lift, lift, @recipe
using LinearAlgebra: normalize

const MaybeObservable{T} = Union{Observable{T}, T}

Makie.plot(f::MaybeObservable{<:AbstractFilament}, args...; kws...) =
    filamentplot(f, args...; kws...)
Makie.plot!(ax, f::MaybeObservable{<:AbstractFilament}, args...; kws...) =
    filamentplot!(ax, f, args...; kws...)

# This macro defines the functions `filamentplot` and `filamentplot!` and the type FilamentPlot,
# among other things.
@recipe FilamentPlot (filament::AbstractFilament, velocities::Union{Nothing, AbstractVector}) begin
    "Number of points per filament segment."
    refinement = 1
    marker = :circle
    markersize = 10.0f0
    linewidth = 1.5f0
    linestyle = :solid
    cycle = [:color]  # this gives the default colour cycle (see e.g. docs for Makie.lines)
    color = @inherit color
    colormap = :viridis
    colorrange = Makie.Automatic()
    markercolor = nothing  # by default, this is the same as `color`

    "Plot tangent vectors?"
    tangents = false

    "Plot curvature vectors?"
    curvatures = false

    "Location of vectors within segment (in [0, 1])."
    vectorpos = 0.0

    visible = true

    tangentcolor = nothing
    curvaturecolor = nothing
    velocitycolor = nothing

    periods = (nothing, nothing, nothing)

    "Arrow attributes (see Makie.arrows3d docs)."
    arrows3d = Attributes()
end

Makie.convert_arguments(::Type{<:FilamentPlot}, f::AbstractFilament) = (f, nothing)  # in case we didn't pass a velocities vector (second argument)

Makie.preferred_axis_type(::FilamentPlot) = Makie.Axis3

fixfirst(f::F, x) where {F} = (ys...) -> f(x, ys...)

function Makie.plot!(p::FilamentPlot)
    f = p.filament
    v = (p.velocities === nothing) ? nothing : p.velocities
    Xs_line = @lift _refine_filament_for_plotting($f, $(p.refinement), $(p.periods))
    Makie.lines!(p, p.attributes, Xs_line)
    Xs_nodes = lift(f, p.periods) do f, Ls
        if all(isnothing, Ls)
            nodes(f)
        else
            map(s⃗ -> Filaments.to_main_periodic_cell(s⃗, Ls), nodes(f))
        end
    end
    Makie.scatter!(p, p.attributes, Xs_nodes)
    let
        tangentcolor = @lift something($(p.tangentcolor), $(p.color))
        _plot_tangents!(
            p, f, p.tangents[], p.vectorpos, p.periods, p.arrows3d;
            color = tangentcolor, colormap = p.colormap, colorrange = p.colorrange,
        )
    end
    return nothing
    let
        curvaturecolor = @lift something($p.curvaturecolor, $p.color)
        _plot_curvatures!(
            p, f, p.curvatures[], p.vectorpos, p.periods, p.arrows3d;
            color = curvaturecolor, colormap = p.colormap, colorrange = p.colorrange,
        )
    end
    if v !== nothing
        velocitycolor = @lift something($p.velocitycolor, $p.color)
        _plot_velocities!(
            p, f, v, p.periods, p.arrows3d;
            color = velocitycolor, colormap = p.colormap, colorrange = p.colorrange,
        )
    end
    p
end

function _refine_filament_for_plotting(
        f::ClosedFilament, refinement::Int, Ls::Tuple,
    )
    Xs_nodes = nodes(f)
    refinement ≥ 1 || error("refinement must be ≥ 1")
    # Note: if `Ls` contains values different from `nothing`, then lines may need to
    # be "broken" by inserting NaNs, and we don't know a priori how many points we need to
    # return.
    with_periods = any(!isnothing, Ls)
    T = eltype(Xs_nodes)
    x⃗_nan = fill(NaN, T)::T
    Xs = similar(Xs_nodes, 0)
    N = refinement * length(f) + 1  # the +1 is to close the loop
    sizehint!(Xs, N)  # this length estimate is only true when we don't need to insert NaNs
    subinds = range(0, 1; length = refinement + 1)[1:refinement]
    Lhs = map(L -> L === nothing ? nothing : L / 2, Ls)  # half periods
    x⃗_prev = Filaments.to_main_periodic_cell(first(Xs_nodes), Ls)
    njumps = 0
    for i ∈ eachindex(f)
        for ζ ∈ subinds
            x⃗ = Filaments.to_main_periodic_cell(f(i, ζ), Ls)
            # Check if we "jumped" from the previous position.
            # Insert a NaN point in that case.
            if with_periods
                if Filaments.is_jump(x⃗, x⃗_prev, Lhs)
                    njumps += 1
                    push!(Xs, x⃗_nan)
                end
            end
            push!(Xs, x⃗)
            x⃗_prev = x⃗
        end
    end
    @assert length(Xs) == refinement * length(f) + njumps
    # Close the loop.
    let x⃗ = Filaments.to_main_periodic_cell(f(lastindex(f), 1.0), Ls)
        if !Filaments.is_jump(x⃗, x⃗_prev, Lhs)
            push!(Xs, x⃗)
        end
    end
    Xs
end

function _tangents_for_arrows(f::AbstractFilament, vectorpos::Real, Ls)
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = f(i, vectorpos, UnitTangent())
        for n ∈ eachindex(x)
            Xs[n][i] = Filaments.to_main_periodic_cell(x[n], Ls[n])
            Vs[n][i] = v[n]
        end
    end
    (Xs..., Vs...)
end

function _plot_arrows!(p, arrows3d, args...; kwargs...)
    attrs = arrows3d[]::Attributes
    Makie.arrows3d!(p, attrs, args...; kwargs...)
end

function _plot_tangents!(
        p::FilamentPlot, f, tangents::Bool, vectorpos, periods, arrows3d;
        kwargs...,
    )
    tangents || return p
    data = @lift _tangents_for_arrows($f, $vectorpos, $periods)
    args = ntuple(i -> @lift(($data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrows3d, args...; kwargs...)
    p
end

function _curvatures_for_arrows(f::AbstractFilament, vectorpos::Real, Ls)
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = f(i, vectorpos, CurvatureVector())
        # ρ = norm(v)  # curvature (= 1 / R)
        for n ∈ eachindex(x)
            Xs[n][i] = Filaments.to_main_periodic_cell(x[n], Ls[n])
            Vs[n][i] = v[n]
        end
    end
    (Xs..., Vs...)
end

function _plot_curvatures!(
        p::FilamentPlot, f, curvatures::Bool,
        vectorpos, periods, arrows3d;
        kwargs...,
    )
    curvatures || return p
    data = @lift _curvatures_for_arrows($f, $vectorpos, $periods)
    args = ntuple(i -> @lift(($data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrows3d, args...; kwargs...)
    p
end

function _velocities_for_arrows(f::AbstractFilament, v::AbstractVector, Ls)
    length(v) == length(f) || throw(DimensionMismatch("wrong length of velocity vector"))
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f), n ∈ eachindex(Xs)
        Xs[n][i] = Filaments.to_main_periodic_cell(f[i][n], Ls[n])
        Vs[n][i] = v[i][n]
    end
    (Xs..., Vs...)
end

function _plot_velocities!(
        p::FilamentPlot, f, v, periods, arrows3d;
        kwargs...,
    )
    data = @lift _velocities_for_arrows($f, $v, $periods)
    args = ntuple(i -> @lift(($data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrows3d, args...; kwargs...)
    p
end

end
