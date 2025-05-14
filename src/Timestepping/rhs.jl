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
    scheduler = DynamicScheduler()  # for threading
    # (0) Apply affect_t! function (optional)
    iter.affect_t!(iter, t)
    # (1) Compute local tangents (this may be used by the forcing)
    _update_unit_tangents!(iter.tangents, fs; scheduler)
    # (2) Compute velocities and optionally streamfunction values
    _update_values_at_nodes!(component, iter.fast_term, fields, fs, t, iter)
    # (3) Modify final velocities if running in non-default mode (optional)
    if iter.mode === MinimalEnergy() && haskey(fields, :velocity)
        # Replace velocities (which will be used for advection) with -s⃗′ × v⃗ₛ
        _minimal_energy_velocities!(fields.velocity, iter.tangents; scheduler)
    end
    fields
end

# Case where only the velocity is passed (generally used in internal RK substeps).
function update_values_at_nodes!(vL::VectorOfVectors, fs, t::Real, iter; kws...)
    fields = (velocity = vL,)
    update_values_at_nodes!(fields, fs, t, iter; kws...)
    vL
end

function _update_unit_tangents!(tangents_all, fs; scheduler)
    for n in eachindex(fs, tangents_all)
        f = fs[n]
        tangents = tangents_all[n]
        tforeach(eachindex(f, tangents); scheduler) do i
            tangents[i] = f[i, UnitTangent()]
        end
    end
    tangents_all
end

function _minimal_energy_velocities!(vs_all, tangents_all; scheduler)
    for n in eachindex(vs_all, tangents_all)
        vs = vs_all[n]
        tangents = tangents_all[n]
        tforeach(eachindex(vs, tangents); scheduler) do i
            vs[i] = vs[i] × tangents[i]  # = -s⃗′ × v⃗ₛ
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
    add_external_fields!(fields, iter, fs, t, iter.to)
    apply_forcing!(fields, iter, fs, t, iter.to)
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
    add_external_fields!(fields, iter, fs, t, iter.to)
    apply_forcing!(fields, iter, fs, t, iter.to)
    nothing
end
