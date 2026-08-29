export Hasimoto

using FFTW
using StaticArrays: SMatrix, SVector
using LinearAlgebra
using .Filaments

struct Hasimoto{Order} <: TemporalScheme end
@inline Hasimoto(order::Int) = Hasimoto{order}()
Hasimoto() = Hasimoto{4}()  # order 4 by default

get_order(::Hasimoto{Order}) where {Order} = Order::Int

function _check_nsubsteps(::Hasimoto, nsubsteps)
    nsubsteps == 1 || @warn("Splittings scheme: the `nsubsteps` parameter is ignored with Hasimoto-based timestepping for local part")
    nothing
end

nbuf_filaments(::Hasimoto) = 0
nbuf_velocities(::Hasimoto) = 0

function s_derivatives(f, Lη, Nf, ks)
    ηs = range(0, Lη; length = Nf+1)
    # Calcul du ŝ
    s_equi = zeros(Nf, 3)
    for i in 1:Nf
        s_equi[i, :] = f(ηs[i])
    end
    s_hat = fft(s_equi, 1) / Nf 
    s_hat[Nf ÷ 2 + 1, :] .= 0

    # Calcul de ŝ', ŝ" et ŝ"'
    s_prime_hat = im .* ks .* s_hat
    s_sec_hat = im .* ks .* s_prime_hat
    s_ter_hat = im .* ks .* s_sec_hat

    # Calcul de s', s" et s"'
    s_prime = real(bfft(s_prime_hat, 1))
    s_sec = real(bfft(s_sec_hat, 1))
    s_ter = real(bfft(s_ter_hat, 1))

    return s_prime, s_sec, s_ter, ηs
end

function curvature_torsion(s_prime, s_sec, s_ter, Nf)
    # Calcul de ρ, τ
    ρ = zeros(Float64, Nf)
    τ = zeros(Float64, Nf)
    for i in 1:Nf
        sp, spp, st = s_prime[i,:], s_sec[i,:], s_ter[i,:]
        cp = sp × spp
        ρ[i] = norm(cp)
        τ[i] = dot(cp, st) / norm(cp)^2
    end
    return ρ, τ
end

function hasimoto_function(ρ, τ, ηs, ks, Nf)
    # Calcul de θ
    func_hat = fft(τ) / Nf
    func_hat[end ÷ 2 + 1] = 0
    θ_hat = @. func_hat / (1im * ks)
    θ_hat[1] = 0.0
    moy = real(func_hat[1]) 
    θ_per = real(bfft(θ_hat))
    θ_per .-= θ_per[1]
    θ = @. θ_per + ηs[1:Nf] * moy

    # θ = similar(τ)
    # @assert length(ηs) == length(τ) + 1
    # L = ηs[end]
    # θ[1] = τ[1] * (ηs[2] - ηs[end] + L) / 2
    # for i in eachindex(τ)[2:end]
    #     θ[i] = θ[i - 1] + τ[i] * (ηs[i + 1] - ηs[i - 1]) / 2
    # end
    # τ_integral = θ[end]
    # τ_mean = τ_integral / L
    # @. θ = θ - θ[1]
    # θ_per = @. θ - ηs[1:Nf] * τ_mean
    # moy = τ_mean

    # Calcul de ψ
    ψ_per = @. ρ * cis(θ_per)
    return ψ_per, θ, moy
end

function orthonormal_frame(s_prime, θ, ks, Nf)
    spp_fourier = fft(s_prime, 1) .* (1im .* ks ./ Nf)
    spp_fourier[Nf ÷ 2 + 1, :] .= 0
    spp = real(bfft(spp_fourier, 1))

    t_hat = zeros(Nf, 3)
    e1_hat = zeros(Nf, 3)
    e2_hat = zeros(Nf, 3)

    @views for i in axes(s_prime, 1)
        s′ = Vec3(s_prime[i, :])
        s″ = Vec3(spp[i, :])
        t̂ = s′ / norm(s′)
        n̂ = s′ × (s″ × s′)
        n̂ = n̂ / norm(n̂)
        b̂ = t̂ × n̂
        t_hat[i, :] = t̂
        e1_hat[i, :] = @. cos(θ[i]) * n̂ - sin(θ[i]) * b̂
        e2_hat[i, :] = @. sin(θ[i]) * n̂ + cos(θ[i]) * b̂
    end

    return t_hat, e1_hat, e2_hat
end

function construct_psi_and_frame(f, Lη, Nf, ks)
    s_prime, s_sec, s_ter, ηs = s_derivatives(f, Lη, Nf, ks)
    s0 = f[1]
    ρ, τ = curvature_torsion(s_prime, s_sec, s_ter, Nf)
    ψper, θ, moy = hasimoto_function(ρ, τ, ηs, ks, Nf)
    t_hat, e1_hat, e2_hat = orthonormal_frame(s_prime, θ, ks, Nf)
    return ψper, moy, t_hat, e1_hat, e2_hat, s0, ηs
end

function psi_and_derivative(ϕ_hat, c, t, k)
    ψper_hat = @. cis(-(k + c)^2 * t) * ϕ_hat
    ψper = bfft(ψper_hat)
    ψp_hat = @. 1im * (k + c) * ψper_hat
    ψp_per = bfft(ψp_hat)
    return ψper, ψp_per
end

######################################

function nls_fourier_nonlinear_if(ϕ_hat, c, t, k)
    ψper, ψp_per = psi_and_derivative(ϕ_hat, c, t, k)
    mb_non_lin = @. 1im / 2 * abs2(ψper) * ψper
    mb_non_lin_hat = fft(mb_non_lin) / length(mb_non_lin)
    # return @. cis((k + c)^2 * t) * mb_non_lin_hat  # XXX: this is not needed, right?
    return mb_non_lin_hat
end

function advance_frame(Nf, T, e1, e2, ψ_per, ψp_per, moy, ηs)
    dT = zeros(Nf, 3)
    de1 = zeros(Nf, 3)
    de2 = zeros(Nf, 3)
    ψ = @. ψ_per * cis(moy * ηs[1:Nf])
    ψp = @. ψp_per * cis(moy * ηs[1:Nf])
    for i in 1:Nf
        a, b, c = real(ψp[i]), imag(ψp[i]), abs2(ψ[i]) / 2
        A = [0.0 -b a ; b 0.0 -c ; -a c 0.0]
        vec = [T[i, :]'; e1[i, :]'; e2[i, :]']
        d_vec = A * vec 

        dT[i, :] = d_vec[1, :]
        de1[i, :] = d_vec[2, :]
        de2[i, :] = d_vec[3, :]
    end
    Tp0 = @. real(ψ[1]) * e1[1, :] + imag(ψ[1]) * e2[1, :]
    ds0 = T[1, :] × Tp0
    return dT, de1, de2, ds0
end

function construct_frame_evolution_matrix(ψ::Complex{T}, ψ′::Complex{T}) where {T}
    a, b, c = real(ψ′), imag(ψ′), abs2(ψ) / 2
    SMatrix{3, 3, T, 9}(
        0, b, -a,  # first *column*
        -b, 0, c,
        a, -c, 0,
    )
end

function RK4IF(ϕ_hat_n, c, T_n, e1_n, e2_n, s0_n, Δt, k, Nf, ηs)
    # Step 1
    tn = zero(Δt)
    ψ1, ψp1 = psi_and_derivative(ϕ_hat_n, c, tn, k)
    v1 = nls_fourier_nonlinear_if(ϕ_hat_n, c, tn, k)
    dT_1, de1_1, de2_1, ds0_1 = advance_frame(Nf, T_n, e1_n, e2_n, ψ1, ψp1, c, ηs)

    # Step 2
    ϕ2 = @. ϕ_hat_n + Δt * v1 / 2
    T_2 = @. T_n .+ Δt .* dT_1 ./ 2
    e1_2 = @. e1_n .+ Δt .* de1_1 ./ 2
    e2_2 = @. e2_n .+ Δt .* de2_1 ./ 2
    t2 = tn + Δt / 2
    ψ2, ψp2 = psi_and_derivative(ϕ2, c, t2, k)
    v2 = nls_fourier_nonlinear_if(ϕ2, c, t2, k)
    dT_2, de1_2, de2_2, ds0_2 = advance_frame(Nf, T_2, e1_2, e2_2, ψ2, ψp2, c, ηs)

    # Step 3
    ϕ3 = @. ϕ_hat_n + Δt * v2 / 2
    T_3 = @. T_n + Δt * dT_2 / 2
    e1_3 = @. e1_n + Δt * de1_2 / 2
    e2_3 = @. e2_n + Δt * de2_2 / 2
    ψ3, ψp3 = psi_and_derivative(ϕ3, c, t2, k)
    v3 = nls_fourier_nonlinear_if(ϕ3, c, t2, k)
    dT_3, de1_3, de2_3, ds0_3 = advance_frame(Nf, T_3, e1_3, e2_3, ψ3, ψp3, c, ηs)

    # Step 4
    ϕ4 = @. ϕ_hat_n + Δt * v3
    T_4 = @. T_n + Δt * dT_3
    e1_4 = @. e1_n + Δt * de1_3
    e2_4 = @. e2_n + Δt * de2_3
    ψ4, ψp4 = psi_and_derivative(ϕ4, c, tn + Δt, k)
    v4 = nls_fourier_nonlinear_if(ϕ4, c, tn + Δt, k)
    dT_4, de1_4, de2_4, ds0_4 = advance_frame(Nf, T_4, e1_4, e2_4, ψ4, ψp4, c, ηs)

    # Final step
    ϕ_hat_np1 = @. ϕ_hat_n + Δt / 6 * (v1 + 2 * v2 + 2 * v3 + v4)
    T_np1 = @. T_n + Δt / 6 * (dT_1 + 2 * dT_2 + 2 * dT_3 + dT_4)
    e1_np1 = @. e1_n + Δt / 6 * (de1_1 + 2 * de1_2 + 2 * de1_3 + de1_4)
    e2_np1 = @. e2_n + Δt / 6 * (de2_1 + 2 * de2_2 + 2 * de2_3 + de2_4)
    s0_np1 = @. s0_n + Δt / 6 * (ds0_1 + 2 * ds0_2 + 2 * ds0_3 + ds0_4)

    return ϕ_hat_np1, T_np1, e1_np1, e2_np1, s0_np1
end

function RK2IF(ϕ_hat_n, c, T_n, e1_n, e2_n, s0_n, Δt, k, Nf, ηs)
    # Step 1
    tn = zero(Δt)
    ψ1, ψp1 = psi_and_derivative(ϕ_hat_n, c, tn, k)
    v1 = nls_fourier_nonlinear_if(ϕ_hat_n, c, tn, k)
    dT_1, de1_1, de2_1, ds0_1 = advance_frame(Nf, T_n, e1_n, e2_n, ψ1, ψp1, c, ηs)

    # Step 2
    ϕ2 = @. ϕ_hat_n + Δt * v1 / 2
    T_2 = @. T_n + Δt * dT_1 / 2
    e1_2 = @. e1_n + Δt * de1_1 / 2
    e2_2 = @. e2_n + Δt * de2_1 / 2
    t2 = tn + Δt / 2
    ψ2, ψp2 = psi_and_derivative(ϕ2, c, t2, k)
    v2 = nls_fourier_nonlinear_if(ϕ2, c, t2, k)
    dT_2, de1_2, de2_2, ds0_2 = advance_frame(Nf, T_2, e1_2, e2_2, ψ2, ψp2, c, ηs)

    # Final step
    ϕ_hat_np1 = @. ϕ_hat_n + Δt * v2
    T_np1 = @. T_n + Δt * dT_2
    e1_np1 = @. e1_n + Δt * de1_2
    e2_np1 = @. e2_n + Δt * de2_2
    s0_np1 = @. s0_n + Δt * ds0_2

    return ϕ_hat_np1, T_np1, e1_np1, e2_np1, s0_np1
end

function NLS_RK2IF(ψ_init_hat, c, Δt, ks)
    ϕ_hat_n = ψ_init_hat

    # Step 1
    tn = zero(Δt)
    v1 = nls_fourier_nonlinear_if(ϕ_hat_n, c, tn, ks)

    # Step 2
    ϕ2 = @. ϕ_hat_n + Δt * v1 / 2
    t2 = tn + Δt / 2
    v2 = nls_fourier_nonlinear_if(ϕ2, c, t2, ks)

    # Final step
    ϕ_hat_np1 = @. ϕ_hat_n + Δt * v2
    ψ_hat_np1 = @. ϕ_hat_np1 * cis(-(ks + c)^2 * Δt)

    return ψ_hat_np1
end

function NLS_Strang2(ψ_init_hat, c, Δt, ks)
    N = length(ψ_init_hat)
    ψ = @. cis(-(ks + c)^2 * Δt/2) * ψ_init_hat  # advance linear term by Δt/2 (Fourier)
    bfft!(ψ)  # to physical space
    @. ψ = cis(Δt * abs2(ψ) / N / 2) * ψ  # advance nonlinear term by Δt (physical space)
    fft!(ψ)  # back to Fourier space
    ψ = @. cis(-(ks + c)^2 * Δt/2) * ψ / N  # advance linear term by Δt/2 (Fourier) + FFT normalisation
    ψ
end

######################################

# Order 2 implementation
function run_hasimoto_simulation(order::Val{2}, f, β, t_in, Δt_in; threshold_ortho = 1e-8)
    Δt = Δt_in * β  # rescale time so we no longer need β
    t = t_in * β    # not sure we need this
    N = length(f)
    Lη = Filaments.knotlims(f)[2]
    # L_actual = integrate(f, Filaments.GaussLegendre(3)) do f, i, ζ
    #     norm(f(i, ζ, Derivative(1)))
    # end
    # @show (Lη - L_actual) / L_actual
    Nf = nextpow(2, N) * 4
    ks = fftfreq(Nf, 2 * π * Nf / Lη)
    ψ_init, moy, T_init, e1_init, e2_init, s0_init, ηs = construct_psi_and_frame(f, Lη, Nf, ks)
    ψ_init_hat = fft(ψ_init) / Nf

    T, e1, e2, s0 = copy(T_init), copy(e1_init), copy(e2_init), copy(s0_init)

    # 1. Advance NLS: ψ(0) -> ψ(Δt/2)
    ψ_hat_mid = NLS_RK2IF(ψ_init_hat, moy, Δt/2, ks)
    # ψ_hat_mid_strang = NLS_Strang2(ψ_init_hat, moy, Δt/2, ks)
    # @show norm(ψ_hat_mid - ψ_hat_mid_strang) / norm(ψ_hat_mid_strang)

    ψp_hat_mid = @. im * (ks + moy) * ψ_hat_mid
    ψ_per_mid = bfft(ψ_hat_mid)
    ψp_per_mid = bfft(ψp_hat_mid)
    ψ_mid = @. ψ_per_mid * cis(moy * ηs[1:Nf])
    ψp_mid = @. ψp_per_mid * cis(moy * ηs[1:Nf])

    # 2. Advance orthonormal frame
    for i in eachindex(ψ_mid)
        A = construct_frame_evolution_matrix(ψ_mid[i], ψp_mid[i])
        # A = @SMatrix [0 -b a; b 0 -c; -a c 0]
        Ω = A * Δt  # first term of Magnus expansion (with GL1 quadrature, i.e. evaluation at midpoint)
        R = exp(Ω)  # this is a unitary matrix (since Ω is skew-symmetric)
        t̂_a = SVector{3}(T_init[i, 1:3])
        ê1_a = SVector{3}(e1_init[i, 1:3])
        ê2_a = SVector{3}(e2_init[i, 1:3])
        X_init = SMatrix{3, 3}(t̂_a..., ê1_a..., ê2_a...)'  # each vector is now a _row_ of X
        X = R * X_init  # rotate frame
        T[i, 1:3] .= X[1, 1:3]
        e1[i, 1:3] .= X[2, 1:3]
        e2[i, 1:3] .= X[3, 1:3]
    end

    # 3. Advance point s0 using some kind of midpoint rule.
    let i = 1
        ψ = ψ_mid[i]  # at Δt/2
        t̂_a = SVector{3}(T_init[i, 1:3])
        t̂_b = SVector{3}(T[i, 1:3])
        e1_a = SVector{3}(e1_init[i, 1:3])
        e1_b = SVector{3}(e1[i, 1:3])
        e2_a = SVector{3}(e2_init[i, 1:3])
        e2_b = SVector{3}(e2[i, 1:3])
        t̂ = (t̂_a + t̂_b) / 2  # estimate tangent at Δt/2
        ê₁ = (e1_a + e1_b) / 2
        ê₂ = (e2_a + e2_b) / 2
        ρ⃗ = real(ψ) * ê₁ + imag(ψ) * ê₂
        ds⃗ = Δt * 1 * (t̂ × ρ⃗)
        s0 = s0_init + ds⃗
    end

    Tp = real.(bfft(fft(T, 1) .* (im .* ks ./ Nf), 1))

    s = filament_reconstruction(T, Nf, ks, s0)
    ξ = Filaments.knots(f)
    Δη = Lη / Nf
    s_ξ = zeros(N, 3)
    for j in 1:N
        # i = mod1(floor(Int, ξ[j] / Δη) + 1, Nf)
        i = unsafe_trunc(Int, (ξ[j] / Lη) * Nf) + 1
        ip1 = mod1(i + 1, Nf)
        t_interp = (ξ[j] - ηs[i]) / Δη
        Xs = @views (Vec3(s[i, 1:3]), Vec3(s[ip1, 1:3]))
        Xsp = @views (Vec3(T[i, 1:3]) .* Δη, Vec3(T[ip1, 1:3]) .* Δη)
        Xspp = @views (Vec3(Tp[i, 1:3]) .* Δη^2, Vec3(Tp[ip1, 1:3]) .* Δη^2)
        # s_ξ[j, :] = Filaments.interpolate(HermiteInterpolation{1}(), Derivative{0}(), t_interp, Xs, Xsp)
        s_ξ[j, :] = Filaments.interpolate(HermiteInterpolation{2}(), Derivative{0}(), t_interp, Xs, Xsp, Xspp)
    end

    return s_ξ, N   
end

# Order 4 implementation
function run_hasimoto_simulation(order::Val{4}, f, β, t_in, Δt_in; threshold_ortho = 1e-8)
    Δt = Δt_in * β  # rescale time so we no longer need β
    t = t_in * β    # not sure we need this
    N = length(f)
    Lη = Filaments.knotlims(f)[2]
    Nf = nextpow(2, N) * 4
    ks = fftfreq(Nf, 2 * π * Nf / Lη)
    ψ_init, moy, T_init, e1_init, e2_init, s0_init, ηs = construct_psi_and_frame(f, Lη, Nf, ks)
    ψ_init_hat = fft(ψ_init) / Nf
    ϕ_hat = ψ_init_hat
    T, e1, e2, s0 = copy(T_init), copy(e1_init), copy(e2_init), copy(s0_init)

    ϕ_hat, T, e1, e2, s0 = RK4IF(ϕ_hat, moy, T, e1, e2, s0, Δt, ks, Nf, ηs)
    max_drift = maximum(abs(norm(T[i, :]) - 1) for i in 1:Nf)
    if max_drift > threshold_ortho
        println("Réorthonormalisation")
        for i in 1:Nf
            T[i, :] ./= norm(T[i, :])
            e1[i, :] .-= dot(e1[i, :], T[i, :]) .* T[i, :] 
            e1[i, :] ./= norm(e1[i, :])
            e2[i, :] .= T[i, :] × e1[i, :]
        end
    end

    ψ_per_hat_final = @. cis(-(ks + moy)^2 * Δt) * ϕ_hat
    ψ_per_final = bfft(ψ_per_hat_final)
    ψ_final = @. ψ_per_final * cis(moy * ηs[1:Nf])

    Tp = real.(bfft(fft(T, 1) .* (im .* ks ./ Nf), 1))

    s = filament_reconstruction(T, Nf, ks, s0)
    ξ = Filaments.knots(f)
    Δη = Lη / Nf
    s_ξ = zeros(N, 3)
    for j in 1:N
        i = unsafe_trunc(Int, (ξ[j] / Lη) * Nf) + 1
        ip1 = mod1(i + 1, Nf)
        t_interp = (ξ[j] - ηs[i]) / Δη
        Xs = @views (Vec3(s[i, 1:3]), Vec3(s[ip1, 1:3]))
        Xsp = @views (Vec3(T[i, 1:3]) .* Δη, Vec3(T[ip1, 1:3]) .* Δη)
        Xspp = @views (Vec3(Tp[i, 1:3]) .* Δη^2, Vec3(Tp[ip1, 1:3]) .* Δη^2)
        # s_ξ[j, :] = Filaments.interpolate(HermiteInterpolation{1}(), Derivative{0}(), t_interp, Xs, Xsp)
        s_ξ[j, :] = Filaments.interpolate(HermiteInterpolation{2}(), Derivative{0}(), t_interp, Xs, Xsp, Xspp)
    end
    return s_ξ, N   
end

function filament_reconstruction(T, Nf, ks, s0)
    Tf = fft(T, 1) / Nf
    sf = @. Tf / (1im * ks)
    sf[1, :] .= 0
    s_rel = real(bfft(sf, 1))
    s = @. s_rel + (s0 - s_rel[1, :])'
    return s
end
