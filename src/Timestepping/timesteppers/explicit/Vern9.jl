# Implementation is adapted from OrdinaryDiffEqVerner.jl (without extra stages or interpolation)
export Vern9

"""
    Vern9 <: ExplicitScheme

Verner's 9-th order, 15 stage Runge–Kutta scheme (Verner, Numer. Algor. *53*, 383, 2010).
"""
struct Vern9 <: ExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::Vern9) = 1
nbuf_velocities(::Vern9) = 14

function _update_velocities!(
        ::Vern9, vs, rhs!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = get_vortex_positions(iter),
    ) where {F <: Function, G <: Function}
    (; fc, vc,) = cache
    fbase = fs
    ftmp = fc[1]
    k1, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15 = vc

    # Note: in principle we don't need to precompute the tableau since all values are
    # compile-time constants.
    tab = Vern9Tableau(number_type(vs), typeof(t))
    (; c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13) = tab
    (;
        a0201, a0301, a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701,
        a0704, a0705, a0706, a0801, a0806, a0807, a0901, a0906, a0907, a0908, a1001, a1006,
        a1007, a1008, a1009, a1101, a1106, a1107, a1108, a1109, a1110, a1201, a1206, a1207,
        a1208, a1209, a1210, a1211, a1301, a1306, a1307, a1308, a1309, a1310, a1311, a1312,
        a1401, a1406, a1407, a1408, a1409, a1410, a1411, a1412, a1413, a1501, a1506, a1507,
        a1508, a1509, a1510, a1511, a1512, a1513, a1514,
        a1601, a1606, a1607, a1608, a1609, a1610, a1611, a1612, a1613,
    ) = tab
    (; b1, b8, b9, b10, b11, b12, b13, b14, b15) = tab

    # Stage 1 (already evaluated)
    @. k1 = vs

    # Stage 2
    # @. vs = a0201 * k1
    advect!(ftmp, vs, a0201 * dt; fbase)
    k2 = vs  # avoid allocation for k2 (we only use it once, in stage 3)
    rhs!(k2, ftmp, t + c1 * dt, iter)

    # Stage 3
    @. vs = a0301 * k1 + a0302 * k2
    advect!(ftmp, vs, dt; fbase)
    rhs!(k3, ftmp, t + c2 * dt, iter)

    # Stage 4
    @. vs = a0401 * k1 + a0403 * k3
    advect!(ftmp, vs, dt; fbase)
    rhs!(k4, ftmp, t + c3 * dt, iter)

    # Stage 5
    @. vs = a0501 * k1 + a0503 * k3 + a0504 * k4
    advect!(ftmp, vs, dt; fbase)
    rhs!(k5, ftmp, t + c4 * dt, iter)

    # Stage 6
    @. vs = a0601 * k1 + a0604 * k4 + a0605 * k5
    advect!(ftmp, vs, dt; fbase)
    rhs!(k6, ftmp, t + c5 * dt, iter)

    # Stage 7
    @. vs = a0701 * k1 + a0704 * k4 + a0705 * k5 + a0706 * k6
    advect!(ftmp, vs, dt; fbase)
    rhs!(k7, ftmp, t + c6 * dt, iter)

    # Stage 8
    @. vs = a0801 * k1 + a0806 * k6 + a0807 * k7
    advect!(ftmp, vs, dt; fbase)
    rhs!(k8, ftmp, t + c7 * dt, iter)

    # Stage 9
    @. vs = a0901 * k1 + a0906 * k6 + a0907 * k7 + a0908 * k8
    advect!(ftmp, vs, dt; fbase)
    rhs!(k9, ftmp, t + c8 * dt, iter)

    # Stage 10
    @. vs = a1001 * k1 + a1006 * k6 + a1007 * k7 + a1008 * k8 + a1009 * k9
    advect!(ftmp, vs, dt; fbase)
    rhs!(k10, ftmp, t + c9 * dt, iter)

    # Stage 11
    @. vs = a1101 * k1 + a1106 * k6 + a1107 * k7 + a1108 * k8 + a1109 * k9 + a1110 * k10
    advect!(ftmp, vs, dt; fbase)
    rhs!(k11, ftmp, t + c10 * dt, iter)

    # Stage 12
    @. vs = a1201 * k1 + a1206 * k6 + a1207 * k7 + a1208 * k8 + a1209 * k9 + a1210 * k10 +
        a1211 * k11
    advect!(ftmp, vs, dt; fbase)
    rhs!(k12, ftmp, t + c11 * dt, iter)

    # Stage 13
    @. vs = a1301 * k1 + a1306 * k6 + a1307 * k7 + a1308 * k8 + a1309 * k9 + a1310 * k10 +
        a1311 * k11 + a1312 * k12
    advect!(ftmp, vs, dt; fbase)
    rhs!(k13, ftmp, t + c12 * dt, iter)

    # Stage 14
    @. vs = a1401 * k1 + a1406 * k6 + a1407 * k7 + a1408 * k8 + a1409 * k9 + a1410 * k10 +
        a1411 * k11 + a1412 * k12 + a1413 * k13
    advect!(ftmp, vs, dt; fbase)
    rhs!(k14, ftmp, t + c13 * dt, iter)

    # Stage 15
    @. vs = a1501 * k1 + a1506 * k6 + a1507 * k7 + a1508 * k8 + a1509 * k9 + a1510 * k10 +
        a1511 * k11 + a1512 * k12 + a1513 * k13 + a1514 * k14
    advect!(ftmp, vs, dt; fbase)
    c14 = 1
    rhs!(k15, ftmp, t + c14 * dt, iter)

    # Stage 16 (for error estimation)
    # @. vs = a1601 * k1 + a1606 * k6 + a1607 * k7 + a1608 * k8 + a1609 * k9 + a1610 * k10 +
    #     a1611 * k11 + a1612 * k12 + a1613 * k13
    # advect!(ftmp, vs, dt; fbase)
    # c15 = 1
    # rhs!(k16, ftmp, t + c15 * dt, iter)

    # Final velocity
    @. vs = b1 * k1 + b8 * k8 + b9 * k9 + b10 * k10 + b11 * k11 + b12 * k12 + b13 * k13 +
        b14 * k14 + b15 * k15

    vs
end

struct Vern9Tableau{T, T2}
    c1::T2
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
    c12::T2
    c13::T2
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
    a0806::T
    a0807::T
    a0901::T
    a0906::T
    a0907::T
    a0908::T
    a1001::T
    a1006::T
    a1007::T
    a1008::T
    a1009::T
    a1101::T
    a1106::T
    a1107::T
    a1108::T
    a1109::T
    a1110::T
    a1201::T
    a1206::T
    a1207::T
    a1208::T
    a1209::T
    a1210::T
    a1211::T
    a1301::T
    a1306::T
    a1307::T
    a1308::T
    a1309::T
    a1310::T
    a1311::T
    a1312::T
    a1401::T
    a1406::T
    a1407::T
    a1408::T
    a1409::T
    a1410::T
    a1411::T
    a1412::T
    a1413::T
    a1501::T
    a1506::T
    a1507::T
    a1508::T
    a1509::T
    a1510::T
    a1511::T
    a1512::T
    a1513::T
    a1514::T
    a1601::T
    a1606::T
    a1607::T
    a1608::T
    a1609::T
    a1610::T
    a1611::T
    a1612::T
    a1613::T
    b1::T
    b8::T
    b9::T
    b10::T
    b11::T
    b12::T
    b13::T
    b14::T
    b15::T
    # btilde1::T
    # btilde8::T
    # btilde9::T
    # btilde10::T
    # btilde11::T
    # btilde12::T
    # btilde13::T
    # btilde14::T
    # btilde15::T
    # btilde16::T
end

Base.@assume_effects :foldable function Vern9Tableau(::Type{T}, ::Type{T2}) where {T <: AbstractFloat, T2 <: AbstractFloat}
    c1 = convert(T2, 0.03462)
    c2 = convert(T2, 0.09702435063878045)
    c3 = convert(T2, 0.14553652595817068)
    c4 = convert(T2, 0.561)
    c5 = convert(T2, 0.22900791159048503)
    c6 = convert(T2, 0.544992088409515)
    c7 = convert(T2, 0.645)
    c8 = convert(T2, 0.48375)
    c9 = convert(T2, 0.06757)
    c10 = convert(T2, 0.25)
    c11 = convert(T2, 0.6590650618730999)
    c12 = convert(T2, 0.8206)
    c13 = convert(T2, 0.9012)
    a0201 = convert(T, 0.03462)
    a0301 = convert(T, -0.03893354388572875)
    a0302 = convert(T, 0.13595789452450918)
    a0401 = convert(T, 0.03638413148954267)
    a0403 = convert(T, 0.10915239446862801)
    a0501 = convert(T, 2.0257639143939694)
    a0503 = convert(T, -7.638023836496291)
    a0504 = convert(T, 6.173259922102322)
    a0601 = convert(T, 0.05112275589406061)
    a0604 = convert(T, 0.17708237945550218)
    a0605 = convert(T, 0.0008027762409222536)
    a0701 = convert(T, 0.13160063579752163)
    a0704 = convert(T, -0.2957276252669636)
    a0705 = convert(T, 0.08781378035642955)
    a0706 = convert(T, 0.6213052975225274)
    a0801 = convert(T, 0.07166666666666667)
    a0806 = convert(T, 0.33055335789153195)
    a0807 = convert(T, 0.2427799754418014)
    a0901 = convert(T, 0.071806640625)
    a0906 = convert(T, 0.3294380283228177)
    a0907 = convert(T, 0.1165190029271823)
    a0908 = convert(T, -0.034013671875)
    a1001 = convert(T, 0.04836757646340646)
    a1006 = convert(T, 0.03928989925676164)
    a1007 = convert(T, 0.10547409458903446)
    a1008 = convert(T, -0.021438652846483126)
    a1009 = convert(T, -0.10412291746271944)
    a1101 = convert(T, -0.026645614872014785)
    a1106 = convert(T, 0.03333333333333333)
    a1107 = convert(T, -0.1631072244872467)
    a1108 = convert(T, 0.03396081684127761)
    a1109 = convert(T, 0.1572319413814626)
    a1110 = convert(T, 0.21522674780318796)
    a1201 = convert(T, 0.03689009248708622)
    a1206 = convert(T, -0.1465181576725543)
    a1207 = convert(T, 0.2242577768172024)
    a1208 = convert(T, 0.02294405717066073)
    a1209 = convert(T, -0.0035850052905728597)
    a1210 = convert(T, 0.08669223316444385)
    a1211 = convert(T, 0.43838406519683376)
    a1301 = convert(T, -0.4866012215113341)
    a1306 = convert(T, -6.304602650282853)
    a1307 = convert(T, -0.2812456182894729)
    a1308 = convert(T, -2.679019236219849)
    a1309 = convert(T, 0.5188156639241577)
    a1310 = convert(T, 1.3653531876033418)
    a1311 = convert(T, 5.8850910885039465)
    a1312 = convert(T, 2.8028087862720628)
    a1401 = convert(T, 0.4185367457753472)
    a1406 = convert(T, 6.724547581906459)
    a1407 = convert(T, -0.42544428016461133)
    a1408 = convert(T, 3.3432791530012653)
    a1409 = convert(T, 0.6170816631175374)
    a1410 = convert(T, -0.9299661239399329)
    a1411 = convert(T, -6.099948804751011)
    a1412 = convert(T, -3.002206187889399)
    a1413 = convert(T, 0.2553202529443446)
    a1501 = convert(T, -0.7793740861228848)
    a1506 = convert(T, -13.937342538107776)
    a1507 = convert(T, 1.2520488533793563)
    a1508 = convert(T, -14.691500408016868)
    a1509 = convert(T, -0.494705058533141)
    a1510 = convert(T, 2.2429749091462368)
    a1511 = convert(T, 13.367893803828643)
    a1512 = convert(T, 14.396650486650687)
    a1513 = convert(T, -0.79758133317768)
    a1514 = convert(T, 0.4409353709534278)
    a1601 = convert(T, 2.0580513374668867)
    a1606 = convert(T, 22.357937727968032)
    a1607 = convert(T, 0.9094981099755646)
    a1608 = convert(T, 35.89110098240264)
    a1609 = convert(T, -3.442515027624454)
    a1610 = convert(T, -4.865481358036369)
    a1611 = convert(T, -18.909803813543427)
    a1612 = convert(T, -34.26354448030452)
    a1613 = convert(T, 1.2647565216956427)
    b1 = convert(T, 0.014611976858423152)
    b8 = convert(T, -0.3915211862331339)
    b9 = convert(T, 0.23109325002895065)
    b10 = convert(T, 0.12747667699928525)
    b11 = convert(T, 0.2246434176204158)
    b12 = convert(T, 0.5684352689748513)
    b13 = convert(T, 0.058258715572158275)
    b14 = convert(T, 0.13643174034822156)
    b15 = convert(T, 0.030570139830827976)
    # bhat1  =convert(T,0.01996996514886773)
    # bhat8  =convert(T,2.19149930494933)
    # bhat9  =convert(T,0.08857071848208438)
    # bhat10 =convert(T,0.11405602348659656)
    # bhat11 =convert(T,0.2533163805345107)
    # bhat12 =convert(T,-2.056564386240941)
    # bhat13 =convert(T,0.340809679901312)
    # bhat16 =convert(T,0.04834231373823958)
    # btilde1 = convert(T, -0.005357988290444578)
    # btilde8 = convert(T, -2.583020491182464)
    # btilde9 = convert(T, 0.14252253154686625)
    # btilde10 = convert(T, 0.013420653512688676)
    # btilde11 = convert(T, -0.02867296291409493)
    # btilde12 = convert(T, 2.624999655215792)
    # btilde13 = convert(T, -0.2825509643291537)
    # btilde14 = convert(T, 0.13643174034822156)
    # btilde15 = convert(T, 0.030570139830827976)
    # btilde16 = convert(T, -0.04834231373823958)

    Vern9Tableau(
        c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, a0201, a0301,
        a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701,
        a0704, a0705, a0706, a0801, a0806, a0807, a0901, a0906, a0907, a0908,
        a1001, a1006, a1007, a1008, a1009, a1101, a1106, a1107, a1108, a1109,
        a1110, a1201, a1206, a1207, a1208, a1209, a1210, a1211, a1301, a1306,
        a1307, a1308, a1309, a1310, a1311, a1312, a1401, a1406, a1407, a1408,
        a1409, a1410, a1411, a1412, a1413, a1501, a1506, a1507, a1508, a1509,
        a1510, a1511, a1512, a1513, a1514, a1601, a1606, a1607, a1608, a1609,
        a1610, a1611, a1612, a1613, b1, b8, b9, b10, b11, b12, b13, b14, b15,
        # btilde1, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13,
        # btilde14, btilde15, btilde16
    )
end
