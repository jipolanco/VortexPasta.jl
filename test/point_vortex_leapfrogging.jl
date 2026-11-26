using Test
using VortexPasta
using VortexPasta.BiotSavart
using VortexPasta.BiotSavart: Vec2
using VortexPasta.Timestepping
using VortexPasta.Diagnostics

function define_leapfrogging_point_vortices(dist_x, dist_y, Ls)
    Lx, Ly = Ls
    x₁ = (Lx - dist_x) / 2
    x₂ = (Lx + dist_x) / 2
    y₁ = (Ly - dist_y) / 2
    y₂ = (Ly + dist_y) / 2
    points = [
        Vec2(x₁, y₁), Vec2(x₂, y₁),  # negative
        Vec2(x₁, y₂), Vec2(x₂, y₂),  # positive
    ]
    signs = [DOWN, DOWN, UP, UP]
    PointVortexState(points, signs)
end

##

T = Float64
Lx, Ly = T(4π), T(2π)
Ls = (Lx, Ly)
Γ = 2.0

backend = CPU()
backend_short = NaiveShortRangeBackend()
backend_long = NonuniformFFTsBackend(backend; σ = 2.0, m = HalfSupport(4))
splitting = @inferred GaussianSplitting(; Ls, rcut = min(Ls...) / 2, β = 4.0)
params = @inferred ParamsPointVortex(; Γ, splitting, backend_short, backend_long)

dist_x = Lx / 5
dist_y = Ly / 5
vortices = define_leapfrogging_point_vortices(dist_x, dist_y, Ls)

##

# Note: for some reason the system destabilises (loses its symmetry) after t ≈ 150.
# This doesn't seem to depend much on the splitting precision, NUFFT or the timestep.
Tsim = 100.0
prob = @inferred PointVortexProblem(vortices, Tsim, params)
iter = @inferred init(prob, Vern8(); dt = 0.4)  # note: this allows _much_ larger timesteps than RK4 to get energy conservation
# iter = @inferred init(prob, RK4(); dt = 0.001)
E₀ = @inferred Diagnostics.kinetic_energy(iter)
p⃗₀ = @inferred Diagnostics.vortex_impulse(iter)
energy = [E₀]
impulse = [p⃗₀]

# using GLMakie
# points_obs = Observable(iter.vortices.points)
# ψs_obs = Observable(iter.ψs)
# time_obs = Observable(iter.t)
# E_obs = Observable(E₀)
#
# fig = Figure()
# ax = Axis(fig[1, 1]; aspect = DataAspect(), title = lift(t -> "t = $t", time_obs))
# text!(ax, 0, Ly; text = lift(E -> "E = $E", E_obs), align = (:left, :top), offset = (2, -2))
# hidespines!(ax)
# hidedecorations!(ax)
# wireframe!(ax, Rect(0, 0, prob.p.Ls...); color = :gray)
# signs = BiotSavart.to_real.(iter.vortices.signs)
# sc = scatter!(ax, points_obs; color = signs, markersize = 8, colormap = :RdBu_4)
# Colorbar(fig[2, 1], sc; vertical = false)
# fig

while iter.t < iter.prob.tspan[2]
    step!(iter)
    if iter.nstep % 10 == 0
        E = Diagnostics.kinetic_energy(iter)
        p⃗ = Diagnostics.vortex_impulse(iter)
        push!(energy, E)
        push!(impulse, p⃗)
        # time_obs[] = iter.t
        # E_obs[] = E
        # notify(ψs_obs)
        # notify(points_obs)  # update plot
        # @show iter.nstep, iter.t, iter.dt, (E - E₀) / E₀
    end
end

@testset "Energy and impulse conservation" begin
    # @show (energy[end] - energy[begin]) / energy[begin]
    @test energy[end] ≈ energy[begin] rtol=1e-6
    @test impulse[end] ≈ impulse[begin] rtol=1e-13
end

## Test annihilation when a dipole gets too close

@testset "Annihilation" begin
    annihilate = AnnihilateFast(0.5 * dist_y)
    iter = @inferred init(prob, Vern8(); dt = 0.4, annihilate)
    @assert length(iter.vortices) == 4
    @assert sum(iter.vortices.signs) == 0

    while iter.t < iter.prob.tspan[2]
        # @show iter.nstep, nvortices, Diagnostics.kinetic_energy(iter)
        step!(iter)
        if length(iter.vortices) < 4
            break  # annihilation happened
        end
    end

    @test length(iter.vortices) == 2     # a dipole was removed
    @test sum(iter.vortices.signs) == 0  # we still have neutral charge
    @test iter.vs[1] ≈ iter.vs[2]  # the two vortices left have the same velocity (they are a dipole, they move in parallel)
end

## Test checkpoints

@testset "Checkpoints" begin
    tsim = 10.0
    T = eltype(vortices)
    prob1 = PointVortexProblem(vortices, tsim, params)
    iter1 = init(prob1, Vern8(); dt = 0.4)
    solve!(iter1)
    checkpoint = mktemp() do path, io
        save_checkpoint(path, iter1)
        checkpoint = @inferred load_checkpoint(path, T)
        @test checkpoint.vortices == iter1.vortices
        @test checkpoint.time == iter1.time
        checkpoint
    end
    @testset "Writing and reading fields" begin
        vs_read = similar(iter1.vs)
        ψs_read = similar(iter1.ψs)
        mktemp() do path, _
            save_checkpoint(path, iter1) do io
                io["velocity"] = iter1.vs
                io["streamfunction"] = iter1.ψs
            end
            local checkpoint = load_checkpoint(path, T) do io
                read!(io, vs_read, "velocity")
                read!(io, ψs_read, "streamfunction")
            end
            @test vs_read == iter1.vs
            @test ψs_read == iter1.ψs
        end
    end
    prob2 = @inferred PointVortexProblem(checkpoint, tsim, params)
    iter2 = @inferred init(prob2, Vern8(); dt = 0.4)
    @test iter2.t == iter1.t ≈ prob1.tspan[2]
    @test iter2.nstep == iter1.nstep > 0
    @test iter2.vortices == iter1.vortices
    @test iter2.vs ≈ iter1.vs  # same velocity (recomputed)
    @test iter2.ψs ≈ iter1.ψs  # same streamfunction (recomputed)
end
