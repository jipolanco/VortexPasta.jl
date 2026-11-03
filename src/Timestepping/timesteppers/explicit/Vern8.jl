# Implementation is adapted from OrdinaryDiffEqVerner.jl (without extra stages or interpolation)
export Vern8

"""
    Vern8 <: ExplicitScheme

Verner's 8-th order, 12 stage Runge–Kutta scheme (Verner, Numer. Algor. *53*, 383, 2010).
"""
struct Vern8 <: ExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::Vern8) = 1
nbuf_velocities(::Vern8) = 11

function _update_velocities!(
        ::Vern8, vs, rhs!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = iter.fs,
    ) where {F <: Function, G <: Function}
    (; fc, vc,) = cache
    fbase = fs
    ftmp = fc[1]
    k1, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12 = vc

    # Note: in principle we don't need to precompute the tableau since all values are
    # compile-time constants.
    tab = Vern8Tableau(number_type(vs), typeof(t))
    (; c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = tab
    (;
        a0201, a0301, a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701,
        a0704, a0705, a0706, a0801, a0804, a0805, a0806, a0807, a0901, a0904, a0905, a0906,
        a0907, a0908, a1001, a1004, a1005, a1006, a1007, a1008, a1009, a1101, a1104, a1105,
        a1106, a1107, a1108, a1109, a1110, a1201, a1204, a1205, a1206, a1207, a1208, a1209,
        a1210, a1211,
        # a1301, a1304, a1305, a1306, a1307, a1308, a1309, a1310,
    ) = tab
    (; b1, b6, b7, b8, b9, b10, b11, b12) = tab

    # Stage 1 (already evaluated)
    @. k1 = vs

    # Stage 2
    # @. vs = a0201 * k1
    advect!(ftmp, vs, a0201 * dt; fbase)
    k2 = vs  # avoid allocation for k2 (we only use it once, in stage 3)
    rhs!(k2, ftmp, t + c2 * dt, iter)

    # Stage 3
    @. vs = a0301 * k1 + a0302 * k2
    advect!(ftmp, vs, dt; fbase)
    rhs!(k3, ftmp, t + c3 * dt, iter)

    # Stage 4
    @. vs = a0401 * k1 + a0403 * k3
    advect!(ftmp, vs, dt; fbase)
    rhs!(k4, ftmp, t + c4 * dt, iter)

    # Stage 5
    @. vs = a0501 * k1 + a0503 * k3 + a0504 * k4
    advect!(ftmp, vs, dt; fbase)
    rhs!(k5, ftmp, t + c5 * dt, iter)

    # Stage 6
    @. vs = a0601 * k1 + a0604 * k4 + a0605 * k5
    advect!(ftmp, vs, dt; fbase)
    rhs!(k6, ftmp, t + c6 * dt, iter)

    # Stage 7
    @. vs = a0701 * k1 + a0704 * k4 + a0705 * k5 + a0706 * k6
    advect!(ftmp, vs, dt; fbase)
    rhs!(k7, ftmp, t + c7 * dt, iter)

    # Stage 8
    @. vs = a0801 * k1 + a0804 * k4 + a0805 * k5 + a0806 * k6 + a0807 * k7
    advect!(ftmp, vs, dt; fbase)
    rhs!(k8, ftmp, t + c8 * dt, iter)

    # Stage 9
    @. vs = a0901 * k1 + a0904 * k4 + a0905 * k5 + a0906 * k6 + a0907 * k7 + a0908 * k8
    advect!(ftmp, vs, dt; fbase)
    rhs!(k9, ftmp, t + c9 * dt, iter)

    # Stage 10
    @. vs = a1001 * k1 + a1004 * k4 + a1005 * k5 + a1006 * k6 + a1007 * k7 + a1008 * k8 +
        a1009 * k9
    advect!(ftmp, vs, dt; fbase)
    rhs!(k10, ftmp, t + c10 * dt, iter)

    # Stage 11
    @. vs = a1101 * k1 + a1104 * k4 + a1105 * k5 + a1106 * k6 + a1107 * k7 + a1108 * k8 +
        a1109 * k9 + a1110 * k10
    advect!(ftmp, vs, dt; fbase)
    rhs!(k11, ftmp, t + c11 * dt, iter)

    # Stage 12
    @. vs = a1201 * k1 + a1204 * k4 + a1205 * k5 + a1206 * k6 + a1207 * k7 + a1208 * k8 +
        a1209 * k9 + a1210 * k10 + a1211 * k11
    advect!(ftmp, vs, dt; fbase)
    c12 = 1
    rhs!(k12, ftmp, t + c12 * dt, iter)

    # Stage 13 (for error estimation)
    # @. vs = a1301 * k1 + a1304 * k4 + a1305 * k5 + a1306 * k6 + a1307 * k7 + a1308 * k8 +
    #     a1309 * k9 + a1310 * k10
    # advect!(ftmp, vs, dt; fbase)
    # c13 = 1
    # rhs!(k13, ftmp, t + c13 * dt, iter)

    # Final velocity
    @. vs = b1 * k1 + b6 * k6 + b7 * k7 + b8 * k8 + b9 * k9 + b10 * k10 + b11 * k11 + b12 * k12

    vs
end

struct Vern8Tableau{T, T2}
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    c7::T2
    c8::T2
    c9::T2
    c10::T2
    c11::T2
    a0201::T
    a0301::T
    a0302::T
    a0401::T
    a0403::T
    a0501::T
    a0503::T
    a0504::T
    a0601::T
    a0604::T
    a0605::T
    a0701::T
    a0704::T
    a0705::T
    a0706::T
    a0801::T
    a0804::T
    a0805::T
    a0806::T
    a0807::T
    a0901::T
    a0904::T
    a0905::T
    a0906::T
    a0907::T
    a0908::T
    a1001::T
    a1004::T
    a1005::T
    a1006::T
    a1007::T
    a1008::T
    a1009::T
    a1101::T
    a1104::T
    a1105::T
    a1106::T
    a1107::T
    a1108::T
    a1109::T
    a1110::T
    a1201::T
    a1204::T
    a1205::T
    a1206::T
    a1207::T
    a1208::T
    a1209::T
    a1210::T
    a1211::T
    a1301::T
    a1304::T
    a1305::T
    a1306::T
    a1307::T
    a1308::T
    a1309::T
    a1310::T
    b1::T
    b6::T
    b7::T
    b8::T
    b9::T
    b10::T
    b11::T
    b12::T
    # btilde1::T
    # btilde6::T
    # btilde7::T
    # btilde8::T
    # btilde9::T
    # btilde10::T
    # btilde11::T
    # btilde12::T
    # btilde13::T
    # extra::Vern8ExtraStages{T, T2}
    # interp::Vern8InterpolationCoefficients{T}
end

Base.@assume_effects :foldable function Vern8Tableau(::Type{T}, ::Type{T2}) where {T <: AbstractFloat, T2 <: AbstractFloat}
    c2 = convert(T2, 0.05)
    c3 = convert(T2, 0.1065625)
    c4 = convert(T2, 0.15984375)
    c5 = convert(T2, 0.39)
    c6 = convert(T2, 0.465)
    c7 = convert(T2, 0.155)
    c8 = convert(T2, 0.943)
    c9 = convert(T2, 0.901802041735857)
    c10 = convert(T2, 0.909)
    c11 = convert(T2, 0.94)
    # c12 = convert(T2, 1)
    # c13 = convert(T2, 1)
    a0201 = convert(T, 0.05)
    a0301 = convert(T, -0.0069931640625)
    a0302 = convert(T, 0.1135556640625)
    a0401 = convert(T, 0.0399609375)
    a0403 = convert(T, 0.1198828125)
    a0501 = convert(T, 0.36139756280045754)
    a0503 = convert(T, -1.3415240667004928)
    a0504 = convert(T, 1.3701265039000352)
    a0601 = convert(T, 0.049047202797202795)
    a0604 = convert(T, 0.23509720422144048)
    a0605 = convert(T, 0.18085559298135673)
    a0701 = convert(T, 0.06169289044289044)
    a0704 = convert(T, 0.11236568314640277)
    a0705 = convert(T, -0.03885046071451367)
    a0706 = convert(T, 0.01979188712522046)
    a0801 = convert(T, -1.767630240222327)
    a0804 = convert(T, -62.5)
    a0805 = convert(T, -6.061889377376669)
    a0806 = convert(T, 5.6508231982227635)
    a0807 = convert(T, 65.62169641937624)
    a0901 = convert(T, -1.1809450665549708)
    a0904 = convert(T, -41.50473441114321)
    a0905 = convert(T, -4.434438319103725)
    a0906 = convert(T, 4.260408188586133)
    a0907 = convert(T, 43.75364022446172)
    a0908 = convert(T, 0.00787142548991231)
    a1001 = convert(T, -1.2814059994414884)
    a1004 = convert(T, -45.047139960139866)
    a1005 = convert(T, -4.731362069449576)
    a1006 = convert(T, 4.514967016593808)
    a1007 = convert(T, 47.44909557172985)
    a1008 = convert(T, 0.01059228297111661)
    a1009 = convert(T, -0.0057468422638446166)
    a1101 = convert(T, -1.7244701342624853)
    a1104 = convert(T, -60.92349008483054)
    a1105 = convert(T, -5.951518376222392)
    a1106 = convert(T, 5.556523730698456)
    a1107 = convert(T, 63.98301198033305)
    a1108 = convert(T, 0.014642028250414961)
    a1109 = convert(T, 0.06460408772358203)
    a1110 = convert(T, -0.0793032316900888)
    a1201 = convert(T, -3.301622667747079)
    a1204 = convert(T, -118.01127235975251)
    a1205 = convert(T, -10.141422388456112)
    a1206 = convert(T, 9.139311332232058)
    a1207 = convert(T, 123.37594282840426)
    a1208 = convert(T, 4.62324437887458)
    a1209 = convert(T, -3.3832777380682018)
    a1210 = convert(T, 4.527592100324618)
    a1211 = convert(T, -5.828495485811623)
    a1301 = convert(T, -3.039515033766309)
    a1304 = convert(T, -109.26086808941763)
    a1305 = convert(T, -9.290642497400293)
    a1306 = convert(T, 8.43050498176491)
    a1307 = convert(T, 114.20100103783314)
    a1308 = convert(T, -0.9637271342145479)
    a1309 = convert(T, -5.0348840888021895)
    a1310 = convert(T, 5.958130824002923)
    b1 = convert(T, 0.04427989419007951)
    b6 = convert(T, 0.3541049391724449)
    b7 = convert(T, 0.24796921549564377)
    b8 = convert(T, -15.694202038838085)
    b9 = convert(T, 25.084064965558564)
    b10 = convert(T, -31.738367786260277)
    b11 = convert(T, 22.938283273988784)
    b12 = convert(T, -0.2361324633071542)
    # bhat1 = convert(T,0.044312615229089795)
    # bhat6 = convert(T,0.35460956423432266)
    # bhat7 = convert(T,0.2478480431366653)
    # bhat8 = convert(T,4.4481347324757845)
    # bhat9 = convert(T,19.846886366118735)
    # bhat10= convert(T,-23.58162337746562)
    # bhat13= convert(T,-0.36016794372897754)
    # btilde1 = convert(T, -3.272103901028138e-5)
    # btilde6 = convert(T, -0.0005046250618777704)
    # btilde7 = convert(T, 0.0001211723589784759)
    # btilde8 = convert(T, -20.142336771313868)
    # btilde9 = convert(T, 5.2371785994398286)
    # btilde10 = convert(T, -8.156744408794658)
    # btilde11 = convert(T, 22.938283273988784)
    # btilde12 = convert(T, -0.2361324633071542)
    # btilde13 = convert(T, 0.36016794372897754)
    #
    # extra = Vern8ExtraStages(T, T2)
    # interp = Vern8InterpolationCoefficients(T)

    Vern8Tableau(
        c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, a0201, a0301, a0302, a0401,
        a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701, a0704, a0705,
        a0706, a0801, a0804, a0805, a0806, a0807, a0901, a0904, a0905, a0906,
        a0907, a0908, a1001, a1004, a1005, a1006, a1007, a1008, a1009, a1101,
        a1104, a1105, a1106, a1107, a1108, a1109, a1110, a1201, a1204, a1205,
        a1206, a1207, a1208, a1209, a1210, a1211, a1301, a1304, a1305, a1306,
        a1307, a1308, a1309, a1310, b1, b6, b7, b8, b9, b10, b11, b12,
        # btilde1, btilde6, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13,
        # extra, interp
    )
end
