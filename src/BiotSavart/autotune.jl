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

!!! warn "Randomness and heuristics"

    The parameters selected by this function can be quite random and change from one run to
    another for a same set of input parameters.
    The heuristics used in this function could be improved and may change in the future.

!!! note "Default accuracy"

    The default values of `β` and `backend_long` are tuned for a nominal 6-digit accuracy.
    See the Extended help for more details.

!!! warn "Periodic domains only"

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

- `ncases = 10`: maximum number of parameter sets to try. Returns the fastest out of all
  tested cases;

- `verbose = false`: if `true`, print autotuning information.

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

_periods_to_tuple(::Type{T}, L::Real) where {T} = ntuple(_ -> convert(T, L), Val(3))
_periods_to_tuple(::Type{T}, L::NTuple{3}) where {T} = map(x -> convert(T, x), L)

function _autotune(
        fs::AbstractVector{<:AbstractFilament{T}}, β::T;
        Ls::NTuple{3, T},
        backend_short = default_short_range_backend(Ls),
        ncases = 10, nruns = 4, verbose = false,
        kws...,
    ) where {T <: AbstractFloat}
    :rcut ∈ keys(kws) && throw(ArgumentError("`rcut` argument must *not* be passed to this function"))
    :Ns ∈ keys(kws) && throw(ArgumentError("`Ns` argument must *not* be passed to this function"))
    :α ∈ keys(kws) && throw(ArgumentError("`α` argument must *not* be passed to this function"))

    Np = sum(length, fs)  # total number of filament nodes
    vs = map(similar ∘ nodes, fs)
    V = prod(Ls)  # domain volume
    α_init = T(1.5) * cbrt(Np / V)  # initial estimate of splitting parameter α
    α::T = α_init

    # Run once and throw away the results to make sure everything is compiled.
    kws_bs = (; Ls, backend_short, kws...,)  # Biot-Savart related kwargs
    _benchmark_params!(vs, fs, α, β; nruns = 1, kws_bs...)

    (; params, score, t, α,) = _benchmark_params!(vs, fs, α, β; nruns, kws_bs...)

    ncase = 0
    verbose && @show ncase, α, t, score

    params_best = params
    times_best = t
    n_best = 0

    for ncase ∈ 1:ncases
        # Select new value of α based on latest score (∈ [-1, 1]), where 0 is ideal.
        # If the score is positive (negative), we need to increase (decrease) α.
        # The score is generally correlated with the actual time spent, but in practice that
        # is not always the case (apparently due to the nontrivial cost of FFTs), so the
        # final decision is done based on times (t < times_best) rather than scores.
        γ = 1 + T(score / 2)  # update factor (heuristic)
        α_prev = α
        α = γ * α
        results = _benchmark_params!(vs, fs, α, β; α_prev, nruns, kws_bs...)  # this is allowed to change α based on backend limits
        results === nothing && break  # can happen if we reached the maximum allowed rcut
        (; params, score, t, α,) = results
        verbose && @show ncase, α, t, score
        if t < times_best
            times_best = t
            params_best = params
            n_best = ncase
        end
        if ncase - n_best ≥ 3
            break  # we're getting too far from the best case, so stop wasting time
        end
    end

    params_best
end

function _benchmark_params!(
        vs, fs, α, β;
        α_prev = zero(α), nruns, Ls, backend_short,
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

    score_best = 1.0
    t_best = Inf

    for _ ∈ 1:nruns
        TimerOutputs.reset_timer!(cache.to)
        t = (@elapsed velocity_on_nodes!(vs, cache, fs))::Float64
        score = time_balance_score(cache)
        if t < t_best
            t_best = t
            score_best = score
        end
    end

    (; params, α, t = t_best, score = score_best,)
end

# Get relative time spent on short and long-range interactions.
#
# Returns a score in [-1, 1], where:
#  * -1 (or 1) mean that too much time is spent on long (or short) range interactions
#  * the value of 0 corresponds to an ideal balance
#
# The implementation depends on the TimerOutputs labels used in the BiotSavart.compute_on_nodes!
# implementations (which differ between CPU-only and CPU+GPU cases).
function time_balance_score(cache::BiotSavartCache)
    device_lr = KA.get_backend(cache.longrange)  # device used for long-range part (CPU or GPU)
    tshort_ratio = _short_range_time_ratio(device_lr, cache)  # in [0, 1]
    2 * tshort_ratio - 1  # in [-1, 1]
end

# CPU-only implementation
# See _compute_on_nodes!(::CPU, ...) for details.
function _short_range_time_ratio(::KA.CPU, cache)
    (; to,) = cache
    tshort = let to = to["Short-range component"]
        # Note: we don't include everything in the "Short-range component" block, but only
        # the most expensive section which scales as N² * rcut³. Other sections don't depend
        # on the parameters we're trying to optimise (α, rcut, kmax).
        TimerOutputs.time(to["Compute Biot–Savart"])
    end
    tlong = let to = to["Long-range component"]
        # Note: as above, we only include the parts which depend on the parameters to
        # optimise. In particular, we don't include the self-interaction correction.
        t = TimerOutputs.time(to["Vorticity to Fourier"]) + TimerOutputs.time(to["Set interpolation points"])
        if haskey(to, "Streamfunction")
            t += TimerOutputs.time(to["Streamfunction"]["Convert to physical"])
        end
        if haskey(to, "Velocity")
            t += TimerOutputs.time(to["Velocity"]["Convert to physical"])
        end
        t
    end
    tshort / (tshort + tlong)
end

# CPU-GPU implementation
# See _compute_on_nodes!(::GPU, ...) for details.
function _short_range_time_ratio(::KA.GPU, cache)
    (; to,) = cache
    (; to_d,) = cache.longrange.common  # timer associated to GPU code running asynchronously
    tshort = let to = to["Short-range component (CPU)"]
        TimerOutputs.time(to["Compute Biot–Savart"])
    end
    tlong_gpu = TimerOutputs.time(to_d["Long-range component (GPU)"])
    tshort / (tshort + tlong_gpu)
end
