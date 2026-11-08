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

# Note: for some reason the system destabilises (loses its symmetry) after t ≈ 150.
# This doesn't seem to depend much on the splitting precision, NUFFT or the timestep.
Tsim = 100.0
prob = @inferred PointVortexProblem(vortices, Tsim, params)
iter = @inferred init(prob, RK4(); dt = 0.04)
E₀ = @inferred Diagnostics.kinetic_energy(iter)
energy = [E₀]

# using GLMakie
# points_obs = Observable(iter.vortices.points)
# ψs_obs = Observable(iter.ψs)
# time_obs = Observable(iter.t)
#
# fig = Figure()
# ax = Axis(fig[1, 1]; aspect = DataAspect(), title = lift(t -> "t = $t", time_obs))
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
        # time_obs[] = iter.t
        # notify(ψs_obs)
        # notify(points_obs)  # update plot
        E = Diagnostics.kinetic_energy(iter)
        push!(energy, E)
        # @show iter.nstep, iter.t, iter.dt, (E - E₀) / E₀
    end
end

@testset "Energy conservation" begin
    @test energy[end] ≈ energy[begin] rtol=1e-6
end
