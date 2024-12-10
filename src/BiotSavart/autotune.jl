# https://discourse.julialang.org/t/is-compat-jl-worth-it-for-the-public-keyword/119041/22
macro public(ex)
    if VERSION >= v"1.11.0-DEV.469"
        args = ex isa Symbol ? (ex,) : Base.isexpr(ex, :tuple) ? ex.args : error("unsupported expression: $ex")
        esc(Expr(:public, args...))
    else
        nothing
    end
end

@public autotune

@doc raw"""
    BiotSavart.autotune(fs::AbstractVector{<:AbstractFilament{T}}, [β::Real = 3.5]; kwargs...) -> ParamsBiotSavart{T}

Generate Biot-Savart parameters optimised for a given filament geometry and accuracy
parameter ``β``.

Based on the given value of ``β``, this function will try to automatically set the following
Ewald-related parameters:

- inverse splitting parameter ``α``;
- cut-off distance ``r_{\text{cut}} = β / α``;
- Fourier-space cut-off ``k_{\text{max}} = 2αβ ≈ πM/L``, where ``M`` is the grid size and
  ``L`` the domain period.

In practice, this function will try to find the value of ``α`` which minimises the
computation of the velocity of the filaments `fs`.

!!! warning "Randomness"

    The parameters selected by this function can be quite random and change from one run to
    another for a same set of input parameters.

!!! note "Default accuracy"

    The default values of `β` and `backend_long` correspond to a nominal 6-digit accuracy.
    See the **Extended help** for more details.

!!! warning "Periodic domains only"

    This function only supports periodic domains (finite domain period ``L``), since the
    chosen parameters are irrelevant in the non-periodic case.

See also [`ParamsBiotSavart`](@ref).

# Extended help

## Mandatory keyword arguments

- `Γ::Real`: vortex circulation (assumed constant);

- `a::Real`: vortex core size (assumed constant);

- `Ls::Union{Real, NTuple{3, Real}}`: domain period in each Cartesian direction.
  If a single value is passed (e.g. `Ls = 2π`), it is assumed that periods are
  the same in all directions.

## Optional keyword arguments

This function accepts the same keyword arguments as [`ParamsBiotSavart`](@ref) and with the
same default values.
In particular, the default short- and long-range backends are:

- `backend_short::ShortRangeBackend = CellListsBackend(2)`;

- `backend_long::LongRangeBackend = NonuniformFFTsBackend(σ = 1.5, m = HalfSupport(4))`.

## Autotuning parameters

The following keyword arguments can be used to control autotuning:

- `nruns = 4`: number of Biot–Savart calculations per value of ``α``. The *minimum* elapsed
  time among all runs will be used in autotuning;

- `Cstart`: initial guess for non-dimensional factor ``C`` (see **Autotuning algorithm** below).
  Default values are 1.5 (pure CPU) and 4.0 (CPU + GPU).

- `ΔC = 0.1`: increment of non-dimensional factor ``C``;

- `verbose = false`: if `true`, print autotuning information.

## Autotuning algorithm

The autotuning algorithm basically consists in trying different values of ``α``, which we
write under the form:

```math
α = C \, {( N / V )}^{1/3}
```

where ``N`` is the total number of filament nodes and ``V`` the domain volume.
The parameter that is varied is the non-dimensional factor ``C``.

For now the algorithm is quite basic.
We try different values around `C = Cstart` using increments of `ΔC`.
The parameters giving the fastest runtime are returned.

## Typical values of ``β`` and NUFFT parameters

The following table roughly relates accuracy (in number of digits) and values of ``β``, as
detailed in [Polanco2024](@citet):

| Precision digits |   ``β``   |  NUFFT ``w`` |
| :--------------: | :-------: | :----------: |
|         3        |    2.0    |      2       |
|         4        |    2.5    |      3       |
|         6        |    3.5    |      4       |
|         8        |    4.0    |      5       |
|        10        |    4.5    |      6       |
|        12        |    5.0    |      7       |
|        14        |    5.5    |      8       |

The last column is the size of the NUFFT half-support ``w`` which ensures sufficient NUFFT accuracy.
The given values assume a NUFFT oversampling factor ``σ = 1.5`` and a (backwards) Kaiser--Bessel
spreading kernel, which are the default when using the [`NonuniformFFTsBackend`](@ref).
Currently, this parameter is not automatically selected by this function.
In other words, knowing the required value of ``w``, one can pass:

    backend_long = NonuniformFFTsBackend(m = HalfSupport(w))

as a keyword argument to this function.

"""
function autotune(
        fs::AbstractVector{<:AbstractFilament{T}},
        β::Real = T(3.5);
        Ls, kws...,
    ) where {T <: AbstractFloat}
    # Convert floats to a common type.
    kws′ = map(values(kws)) do val
        val isa AbstractFloat ? convert(T, val) : val
    end
    Ls′ = _periods_to_tuple(T, Ls)
    _autotune(fs, convert(T, β); kws′..., Ls = Ls′)
end

_periods_to_tuple(::Type{T}, L::Real) where {T} = _periods_to_tuple(T, (L, L, L))
_periods_to_tuple(::Type{T}, L::NTuple{3}) where {T} = map(x -> convert(T, x), L)

# Default value of Cstart: if the long-range part runs on a (fast) GPU, the optimal value of
# α (and C) is expected to be larger than in the pure CPU case.
# TODO: should this depend on the actual GPU device, and in particular on its memory limits?
# (larger α uses more memory)
default_Cstart(::KA.CPU) = 1.5
default_Cstart(::KA.GPU) = 4.0
default_Cstart(backend::LongRangeBackend) = default_Cstart(KA.get_backend(backend))

function _autotune(
        fs::AbstractVector{<:AbstractFilament{T}}, β::T;
        Ls::NTuple{3, T},
        backend_short = default_short_range_backend(Ls),
        backend_long = default_long_range_backend(Ls),
        nruns = 4,
        Cstart::T = T(default_Cstart(backend_long)),
        ΔC::T = T(0.1),
        verbose = false,
        kws...,
    ) where {T <: AbstractFloat}
    :rcut ∈ keys(kws) && throw(ArgumentError("`rcut` argument must *not* be passed to this function"))
    :Ns ∈ keys(kws) && throw(ArgumentError("`Ns` argument must *not* be passed to this function"))
    :α ∈ keys(kws) && throw(ArgumentError("`α` argument must *not* be passed to this function"))

    Np = sum(length, fs)  # total number of filament nodes
    vs = map(similar ∘ nodes, fs)
    V = prod(Ls)  # domain volume

    α_base::T = cbrt(Np / V)  # α = C * α_base

    # Run once and throw away the results to make sure everything is compiled.
    kws_bs = (; Ls, backend_short, backend_long, kws...,)  # Biot-Savart related kwargs
    _benchmark_params!(vs, fs, Cstart * α_base, β; nruns = 1, kws_bs...)

    # Run the base case: C = Cstart.
    C = Cstart
    results_base = _benchmark_params!(vs, fs, C * α_base, β; nruns, kws_bs...)
    @assert results_base !== nothing
    let results = results_base
        verbose && @show C, results.α, results.t
    end
    α_prev = results_base.α
    params_best = results_base.params
    t_best = results_base.t

    # Try larger values of C (and α)
    C = Cstart
    nworse = 0  # how many successive trials are worse than the best one?
    while nworse < 2  # allow 2 successive "worse" cases
        C = C + ΔC
        results = _benchmark_params!(vs, fs, C * α_base, β; α_prev, nruns, kws_bs...)
        # Since we're increasing α and thus reducing rcut, the `continue` won't
        # happen indefinitely, and will stop for some value of α.
        results === nothing && continue
        verbose && @show C, results.α, results.t
        α_prev = results.α
        (; t,) = results
        if t < t_best
            t_best = t
            params_best = results.params
            nworse = 0
        else
            # This case is worse than the best case.
            # But we don't exit immediately, in case the trend is not monotonic.
            nworse += 1
        end
    end

    # Try smaller values of C (and α)
    C = Cstart
    α_prev = results_base.α
    nworse = 0
    while nworse < 2
        C = C - ΔC
        results = _benchmark_params!(vs, fs, C * α_base, β; α_prev, nruns, kws_bs...)
        # Since we're decreasing α and thus increasing rcut, we stop here if we're limited
        # by rcut_max.
        results === nothing && break
        verbose && @show C, results.α, results.t
        α_prev = results.α
        (; t,) = results
        if t < t_best
            t_best = t
            params_best = results.params
            nworse = 0
        else
            nworse += 1
        end
    end

    if verbose
        α_best = params_best.α
        C_best = α_best / α_base
        @info "[BiotSavart.autotune] Best case:" C_best α_best t_best
    end

    params_best
end

function _benchmark_params!(
        vs, fs, α, β;
        α_prev = nothing, nruns, Ls, backend_short,
        timer = TimerOutput(),  # useful for debugging; can be directly passed to autotune
        kws...,
    )
    T = Filaments.number_type(fs)
    @assert T === eltype(Ls)
    rcut_max = max_cutoff_distance(backend_short, Ls)
    rcut = β / α

    # NOTE: we're allowed to adjust α to respect the maximum cut-off distance rcut_max
    # required by the short-range backend.
    if rcut > rcut_max
        rcut = rcut_max
        α = β / rcut
        α == α_prev && return nothing  # don't repeat a previous benchmark if we're limited by rcut_max
    end

    kmax = 2 * α * β
    Ns = map(Ls) do L
        ceil(Int, kmax * L / T(π) + 1)
    end
    params = ParamsBiotSavart(T; Ls, α, rcut, backend_short, Ns, kws...)
    cache = init_cache(params, fs; timer)

    t_best = Inf

    for _ ∈ 1:nruns
        TimerOutputs.reset_timer!(cache.to)
        stats = @timed velocity_on_nodes!(vs, cache, fs)
        # @show stats.gctime, stats.lock_conflicts, stats.compile_time, stats.recompile_time
        t = stats.time::Float64
        if t < t_best
            t_best = t
        end
    end

    (; params, α, t = t_best,)
end
