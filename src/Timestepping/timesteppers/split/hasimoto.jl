export Hasimoto

using FFTW
using LinearAlgebra
using .Filaments

struct Hasimoto <: SplittingScheme end

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
    func = τ
    func_hat = fft(func) / Nf
    func_hat[Nf ÷ 2 + 1] = 0
    θ_hat = @. func_hat / (1im * ks)
    θ_hat[1] = 0.0
    moy = real(func_hat[1]) 
    θ_per = real(bfft(θ_hat))
    θ_per .-= θ_per[1]
    θ = @. θ_per + ηs[1:Nf] * moy

    # Calcul de ψ
    ψ_per = zeros(ComplexF64, Nf)
    @. ψ_per = ρ * cis(θ_per)
    return ψ_per, θ, moy
end

function orthonormal_frame(s_prime, θ, ks, Nf)
    # Calcul de t_hat 
    t_hat = zeros(Nf, 3)
    for i in 1:Nf
        t_hat[i, :] = s_prime[i, :] ./ norm(s_prime[i, :])
    end

    # Calcul de e1_hat et e2_hat
    tf_hat = fft(t_hat, 1) / Nf
    tf_hat[Nf ÷ 2 + 1, :] .= 0
    tfp_hat = @. 1im * ks * tf_hat
    tp_hat = real(bfft(tfp_hat, 1))

    n_hat = zeros(Nf, 3)
    b_hat = zeros(Nf, 3)
    e1_hat = zeros(Nf, 3)
    e2_hat = zeros(Nf, 3)
    for i in 1:Nf
        n_hat[i, :] = tp_hat[i, :] ./ norm(tp_hat[i, :])
        b_hat[i, :] = t_hat[i, :] × n_hat[i, :]
        e1_hat[i, :] = @. cos(θ[i]) * n_hat[i, :] - sin(θ[i]) * b_hat[i, :]
        e2_hat[i, :] = @. sin(θ[i]) * n_hat[i, :] + cos(θ[i]) * b_hat[i, :]
    end

    return t_hat, e1_hat, e2_hat
end

function construct_psi_and_frame(f, Lη, Nf, ks)
    s_prime, s_sec, s_ter, ηs = s_derivatives(f, Lη, Nf, ks)
    s0 = f(ηs[1])
    ρ, τ = curvature_torsion(s_prime, s_sec, s_ter, Nf)
    ψper, θ, moy = hasimoto_function(ρ, τ, ηs, ks, Nf)
    t_hat, e1_hat, e2_hat = orthonormal_frame(s_prime, θ, ks, Nf)
    return ψper, moy, t_hat, e1_hat, e2_hat, s0, ηs
end

function psi_and_derivative(ϕ_hat, c, β, t, k, Nf)
    ψper_hat = @. cis(-β * (k + c) ^ 2 * t) * ϕ_hat
    ψper = bfft(ψper_hat)
    ψp_hat = @. 1im * (k + c) * ψper_hat
    ψp_per = bfft(ψp_hat)
    return ψper, ψp_per
end
######################################

function g(ϕ_hat, c, β, t, k, Nf)
    ψper, ψp_per = psi_and_derivative(ϕ_hat, c, β, t, k, Nf)
    mb_non_lin = @. 1im / 2 * β * abs(ψper) ^ 2 * ψper
    mb_non_lin_hat = fft(mb_non_lin) / Nf
    return @. cis(β * (k + c) ^ 2 * t) * mb_non_lin_hat
end

function h(Nf, T, e1, e2, ψ_per, ψp_per, moy, ηs, β)
    dT = zeros(Nf, 3)
    de1 = zeros(Nf, 3)
    de2 = zeros(Nf, 3)
    ψ = @. ψ_per * cis(moy * ηs[1:Nf])
    ψp = @. ψp_per * cis(moy * ηs[1:Nf])
    for i in 1:Nf
        a, b, c = β * real(ψp[i]), β * imag(ψp[i]), β * abs(ψ[i])^2 / 2
        A = [0.0 -b a ; b 0.0 -c ; -a c 0.0]
        vec = [T[i, :]'; e1[i, :]'; e2[i, :]']
        d_vec = A * vec 

        dT[i, :] = d_vec[1, :]
        de1[i, :] = d_vec[2, :]
        de2[i, :] = d_vec[3, :]
    end
    Tp0 = @. real(ψ[1]) * e1[1, :] + imag(ψ[1]) * e2[1, :]
    ds0 = β .* (T[1, :] × Tp0)
    return dT, de1, de2, ds0
end

function RK4IF(ϕ_hat_n, c, T_n, e1_n, e2_n, s0_n, β, tn, Δt, k, Nf, ηs)
    # Step 1
    ψ1, ψp1 = psi_and_derivative(ϕ_hat_n, c, β, tn, k, Nf)
    v1 = g(ϕ_hat_n, c, β, tn, k, Nf)
    dT_1, de1_1, de2_1, ds0_1 = h(Nf, T_n, e1_n, e2_n, ψ1, ψp1, c, ηs, β)

    # Step 2
    ϕ2 = @. ϕ_hat_n + Δt * v1 / 2
    T_2 = @. T_n .+ Δt .* dT_1 ./ 2
    e1_2 = @. e1_n .+ Δt .* de1_1 ./ 2
    e2_2 = @. e2_n .+ Δt .* de2_1 ./ 2
    t2 = tn + Δt / 2
    ψ2, ψp2 = psi_and_derivative(ϕ2, c, β, t2, k, Nf)
    v2 = g(ϕ2, c, β, t2, k, Nf)
    dT_2, de1_2, de2_2, ds0_2 = h(Nf, T_2, e1_2, e2_2, ψ2, ψp2, c, ηs, β)

    # Step 3
    ϕ3 = @. ϕ_hat_n + Δt * v2 / 2
    T_3 = @. T_n + Δt * dT_2 / 2
    e1_3 = @. e1_n + Δt * de1_2 / 2
    e2_3 = @. e2_n + Δt * de2_2 / 2
    ψ3, ψp3 = psi_and_derivative(ϕ3, c, β, t2, k, Nf)
    v3 = g(ϕ3, c, β, t2, k, Nf)
    dT_3, de1_3, de2_3, ds0_3 = h(Nf, T_3, e1_3, e2_3, ψ3, ψp3, c, ηs, β)

    # Step 4
    ϕ4 = @. ϕ_hat_n + Δt * v3
    T_4 = @. T_n + Δt * dT_3
    e1_4 = @. e1_n + Δt * de1_3
    e2_4 = @. e2_n + Δt * de2_3
    ψ4, ψp4 = psi_and_derivative(ϕ4, c, β, tn + Δt, k, Nf)
    v4 = g(ϕ4, c, β, tn + Δt, k, Nf)
    dT_4, de1_4, de2_4, ds0_4 = h(Nf, T_4, e1_4, e2_4, ψ4, ψp4, c, ηs, β)

    # Final step
    ϕ_hat_np1 = @. ϕ_hat_n + Δt / 6 * (v1 + 2 * v2 + 2 * v3 + v4)
    T_np1 = @. T_n + Δt / 6 * (dT_1 + 2 * dT_2 + 2 * dT_3 + dT_4)
    e1_np1 = @. e1_n + Δt / 6 * (de1_1 + 2 * de1_2 + 2 * de1_3 + de1_4)
    e2_np1 = @. e2_n + Δt / 6 * (de2_1 + 2 * de2_2 + 2 * de2_3 + de2_4)
    s0_np1 = @. s0_n + Δt / 6 * (ds0_1 + 2 * ds0_2 + 2 * ds0_3 + ds0_4)

    return ϕ_hat_np1, T_np1, e1_np1, e2_np1, s0_np1
end

######################################

function run_hasimoto_simulation(f, β, t, Δt; threshold_ortho = 1e-8)
    N = length(f)
    Lη = Filaments.knotlims(f)[2]
    Nf = nextpow(2, N)
    ks = fftfreq(Nf, 2 * π * Nf / Lη)
    ψ_init, moy, T_init, e1_init, e2_init, s0_init, ηs = construct_psi_and_frame(f, Lη, Nf, ks)
    ψ_init_hat = fft(ψ_init) / Nf
    ϕ_hat = ψ_init_hat
    T, e1, e2, s0 = copy(T_init), copy(e1_init), copy(e2_init), copy(s0_init)

    ϕ_hat, T, e1, e2, s0 = RK4IF(ϕ_hat, moy, T, e1, e2, s0, β, t, Δt, ks, Nf, ηs)
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

    ψ_per_hat_final = @. cis(-β * (ks + moy) ^ 2 * (t + Δt)) * ϕ_hat
    ψ_per_final = bfft(ψ_per_hat_final)
    ψ_final = @. ψ_per_final * cis(moy * ηs[1:Nf]) 

    s = filament_reconstruction(T, Nf, ks, s0)
    ξ = Filaments.knots(f)
    Δη = Lη / Nf
    s_ξ = zeros(N, 3)
    for j in 1:N
        i = mod1(floor(Int, ξ[j] / Δη) + 1, Nf)
        ip1 = mod1(i + 1, Nf)
        t_interp = (ξ[j] - ηs[i]) / Δη
        Xs = (s[i, :], s[ip1, :])
        Xsp = (T[i, :] .* Δη, T[ip1, :] .* Δη)
        s_ξ[j, :] = Filaments.interpolate(HermiteInterpolation{1}(), Derivative{0}(), t_interp, Xs, Xsp)
        # println("s_ξ[j] avec interpolation = ", s_ξ[j, :])
        # s_ξ[j, :] = s[j, :]
        # println("s_ξ[j] sans interpolation = ", s_ξ[j, :])
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