# These functions will compute the right-hand side of the evolution equation, i.e. the
# velocity v⃗ in ds⃗/dt = v⃗(s⃗, t). This usually includes the vortex induced velocity
# (Biot–Savart law) as well as possible additional velocities (forcing, ...).

# This is the most general variant which should be called by timesteppers.
# Throughout the Timestepping module, `rhs!` is basically an alias to this function.
#
# Here `fields` can include the quantities:
# - `velocity`: will contain the instantaneous velocity of filament nodes
# - `streamfunction` (optional): if included, will contain streamfunction values at filament nodes
function update_values_at_nodes!(
        fields::NamedTuple, fs, t::Real, iter::VortexFilamentSolver;
        component = Val(:full),  # compute slow + fast components by default
    )
    (; to, quantities) = iter
    # (0) Apply affect_t! function (optional)
    if iter.affect_t! !== default_callback()
        @timeit to "Affect! (fine-grained)" iter.affect_t!(iter, t)
    end
    # (1) Compute local tangents (this may be used by the forcing).
    # In fact, this is currently only used by NormalFluidForcing or in MinimalEnergy mode.
    if hasproperty(quantities, :tangents)
        @timeit to "Compute tangents" _update_unit_tangents!(quantities.tangents, fs; scheduler)
    end
    # (2) Compute velocities and optionally streamfunction values
    @timeit to "Compute filament velocity" _update_values_at_nodes!(component, iter.fast_term, fields, fs, t, iter)
    # (3) Modify final velocities if running in non-default mode (optional)
    if iter.mode === MinimalEnergy() && haskey(fields, :velocity)
        # Replace velocities (which will be used for advection) with -s⃗′ × v⃗ₛ
        @assert hasproperty(quantities, :tangents)
        @timeit to "MinimalEnergy mode" _minimal_energy_velocities!(fields.velocity, quantities.tangents, iter; scheduler)
    end
    fields
end

# Case where only the velocity is passed (generally used in internal RK substeps).
function update_values_at_nodes!(vL::VectorOfVectors, fs, t::Real, iter; kws...)
    fields = (velocity = vL,)
    update_values_at_nodes!(fields, fs, t, iter; kws...)
    vL
end

function _update_unit_tangents!(tangents_all, fs)
    @sync for chunk in FilamentChunkIterator(fs)
        Threads.@spawn for (n, inds, _) in chunk
            f = fs[n]
            tangents = tangents_all[n]
            for i in inds
                @inbounds tangents[i] = f[i, UnitTangent()]
            end
        end
    end
    tangents_all
end

function _minimal_energy_velocities!(vL_all, tangents_all, iter; scheduler)
    vs_all = iter.vs  # BS velocities will be copied here
    for n in eachindex(vL_all, vs_all, tangents_all)
        vL = vL_all[n]
        vs = vs_all[n]
        tangents = tangents_all[n]
        tforeach(eachindex(vL, vs, tangents); scheduler) do i
            vs[i] = vL[i]  # copy actual BS velocity into iter.vs
            vL[i] = vL[i] × tangents[i]  # = -s⃗′ × v⃗ₛ
        end
    end
    vs_all
end

# This variant computes the full BS law + any added velocities.
function _update_values_at_nodes!(
        ::Val{:full},
        ::FastBiotSavartTerm,  # ignored in this case
        fields::NamedTuple{Names, NTuple{N, V}},
        fs::VectorOfFilaments,
        t::Real,
        iter::VortexFilamentSolver,
    ) where {Names, N, V <: VectorOfVectors}
    if iter.LIA
        BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(:only))
    else
        BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs)
    end
    add_external_fields!(fields, iter, fs, t, iter.to)
    apply_forcing!(fields, iter, fs, t, iter.to)
    nothing
end

# Compute slow component only.
# This is generally called in IMEX-RK substeps, where only the velocity (and not the
# streamfunction) is needed.
# We assume that the "slow" component is everything but LIA term when evolving the
# Biot-Savart law.
# This component will be treated explicitly by IMEX schemes.
function _update_values_at_nodes!(
        ::Val{:slow},
        ::LocalTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        t::Real,
        iter::VortexFilamentSolver,
    )
    T = eltype_nested(Vec3, fields.velocity)
    @assert T <: Vec3
    if iter.LIA
        fill!(fields.velocity, zero(T))
    else
        BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(false))
    end
    add_external_fields!(fields, iter, fs, t, iter.to)
    apply_forcing!(fields, iter, fs, t, iter.to)
    nothing
end

function _update_values_at_nodes!(
        ::Val{:slow},
        ::ShortRangeTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        t::Real,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; shortrange = false)
    add_external_fields!(fields, iter, fs, t, iter.to)
    apply_forcing!(fields, iter, fs, t, iter.to)
    nothing
end

# Compute fast component only.
# This is generally called in IMEX-RK substeps, where only the velocity (and not the
# streamfunction) is needed.
# We assume that the "fast" component is the LIA term when evolving the Biot-Savart law.
# This component will be treated implicitly by IMEX schemes.
function _update_values_at_nodes!(
        ::Val{:fast},
        ::LocalTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        t::Real,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(:only))
    nothing
end

function _update_values_at_nodes!(
        ::Val{:fast},
        ::ShortRangeTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        t::Real,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; longrange = false)
    nothing
end
