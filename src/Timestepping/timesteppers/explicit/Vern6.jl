# Implementation is adapted from OrdinaryDiffEqVerner.jl (without extra stages or interpolation)
export Vern6

"""
    Vern6 <: ExplicitScheme

Verner's 6-th order, 9 stage Runge–Kutta scheme (Verner, Numer. Algor. *53*, 383, 2010).
"""
struct Vern6 <: ExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::Vern6) = 1
nbuf_velocities(::Vern6) = 6

function _update_velocities!(
        ::Vern6, vs, rhs!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = iter.fs,
    ) where {F <: Function, G <: Function}
    (; fc, vc,) = cache
    fbase = fs
    ftmp = fc[1]
    k1, k3, k4, k5, k6, k7 = vc

    # Note: in principle we don't need to precompute the tableau since all values are
    # compile-time constants.
    tab = Vern6Tableau(number_type(vs), typeof(t))
    (;
        c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65,
        a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98,
    ) = tab

    # Stage 1 (already evaluated)
    @. k1 = vs

    # Stage 2
    # @. vs = a21 * k1
    advect!(ftmp, vs, a21 * dt; fbase)
    k2 = vs  # avoid allocation for k2 (we only use it once, in stage 3)
    rhs!(k2, ftmp, t + c1 * dt, iter)

    # Stage 3
    @. vs = a31 * k1 + a32 * k2   # we should divide this by c2 to get an actual velocity for propagating from t -> t + c2 * dt...
    advect!(ftmp, vs, dt; fbase)  # but that is equivalent to propagating from t -> t + dt with velocity vs = v_actual * c2 (so we put dt instead of c2 * dt here)
    rhs!(k3, ftmp, t + c2 * dt, iter)

    # Stage 4
    @. vs = a41 * k1 + a43 * k3
    advect!(ftmp, vs, dt; fbase)
    rhs!(k4, ftmp, t + c3 * dt, iter)

    # Stage 5
    @. vs = a51 * k1 + a53 * k3 + a54 * k4
    advect!(ftmp, vs, dt; fbase)
    rhs!(k5, ftmp, t + c4 * dt, iter)

    # Stage 6
    @. vs = a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5
    advect!(ftmp, vs, dt; fbase)
    rhs!(k6, ftmp, t + c5 * dt, iter)

    # Stage 7
    @. vs = a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6
    advect!(ftmp, vs, dt; fbase)
    rhs!(k7, ftmp, t + c6 * dt, iter)

    # Stage 8
    @. vs = a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7
    c7 = 1  # advance to t + dt
    advect!(ftmp, vs, dt; fbase)
    k8 = vs  # avoid allocation for k8 (reuse vs)
    rhs!(k8, ftmp, t + c7 * dt, iter)

    # Stage 9
    @. vs = a91 * k1 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7 + a98 * k8

    vs
end

struct Vern6Tableau{T, T2}
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a43::T
    a51::T
    a53::T
    a54::T
    a61::T
    a63::T
    a64::T
    a65::T
    a71::T
    a73::T
    a74::T
    a75::T
    a76::T
    a81::T
    a83::T
    a84::T
    a85::T
    a86::T
    a87::T
    a91::T
    a94::T
    a95::T
    a96::T
    a97::T
    a98::T
    # btilde1::T
    # btilde4::T
    # btilde5::T
    # btilde6::T
    # btilde7::T
    # btilde8::T
    # btilde9::T
    # extra::Vern6ExtraStages{T, T2}
    # interp::Vern6InterpolationCoefficients{T}
end

function Vern6Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    c1 = convert(T2, 3 // 50)
    c2 = convert(T2, 1439 // 15000)
    c3 = convert(T2, 1439 // 10000)
    c4 = convert(T2, 4973 // 10000)
    c5 = convert(T2, 389 // 400)
    c6 = convert(T2, 1999 // 2000)
    a21 = convert(T, 3 // 50)
    a31 = convert(T, 519479 // 27000000)
    a32 = convert(T, 2070721 // 27000000)
    a41 = convert(T, 1439 // 40000)
    a43 = convert(T, 4317 // 40000)
    a51 = convert(T, 109225017611 // 82828840000)
    a53 = convert(T, -417627820623 // 82828840000)
    a54 = convert(T, 43699198143 // 10353605000)
    a61 = convert(T, -8036815292643907349452552172369 // 191934985946683241245914401600)
    a63 = convert(T, 246134619571490020064824665 // 1543816496655405117602368)
    a64 = convert(T, -13880495956885686234074067279 // 113663489566254201783474344)
    a65 = convert(T, 755005057777788994734129 // 136485922925633667082436)
    a71 = convert(
        T,
        # -BigInt(1663299841566102097180506666498880934230261) //
        #     BigInt(30558424506156170307020957791311384232000)
        -54.430156935316504  # we hardcode the value to make sure it's a compile-time constant (this wouldn't be OK if we used extended precision, BigFloat or similar)
    )
    a73 = convert(
        T,
        130838124195285491799043628811093033 // 631862949514135618861563657970240
    )
    a74 = convert(
        T,
        # -BigInt(3287100453856023634160618787153901962873) //
        #     BigInt(20724314915376755629135711026851409200)
        -158.61081378459
    )
    a75 = convert(
        T,
        2771826790140332140865242520369241 // 396438716042723436917079980147600
    )
    a76 = convert(T, -1799166916139193 // 96743806114007800)
    a81 = convert(
        T,
        # -BigInt(832144750039369683895428386437986853923637763) //
        #     BigInt(15222974550069600748763651844667619945204887)
        -54.66374178728198
    )
    a83 = convert(
        T,
        818622075710363565982285196611368750 //
            3936576237903728151856072395343129
    )
    a84 = convert(
        T,
        # -BigInt(9818985165491658464841194581385463434793741875) //
        #     BigInt(61642597962658994069869370923196463581866011)
        -159.2889574744995
    )
    a85 = convert(
        T,
        # BigInt(31796692141848558720425711042548134769375) //
        #     BigInt(4530254033500045975557858016006308628092)
        7.018743740796944
    )
    a86 = convert(T, -14064542118843830075 // 766928748264306853644)
    a87 = convert(T, -1424670304836288125 // 2782839104764768088217)
    a91 = convert(T, 382735282417 // 11129397249634)
    a94 = convert(T, 5535620703125000 // 21434089949505429)
    a95 = convert(T, 13867056347656250 // 32943296570459319)
    a96 = convert(T, 626271188750 // 142160006043)
    a97 = convert(T, -51160788125000 // 289890548217)
    a98 = convert(T, 163193540017 // 946795234)
    # b1   =convert(T,124310637869885675646798613//2890072468789466426596827670)
    # b4   =convert(T,265863151737164990361330921875//1113197271463372303940319369579)
    # b5   =convert(T,3075493557174030806536302953125//6843749922042323876546949699876)
    # b6   =convert(T,67798000008733879813263055//29532792147666737550036372)
    # b7   =convert(T,-1099436585155390846238326375//15055706496446408859196167)
    # b8   =convert(T,26171252653086373181571802//368794478890732346033505)
    # b9   =convert(T,1//30)
    # btilde1 = convert(T, 12461131651614938103148389 // 1445036234394733213298413835)
    # btilde4 = convert(T, -21633909117387045317965953125 // 1113197271463372303940319369579)
    # btilde5 = convert(T, 21633909117387045317965953125 // 760416658004702652949661077764)
    # btilde6 = convert(T, -6922850917563854501749105 // 3281421349740748616670708)
    # btilde7 = convert(T, 173071272939096362543727625 // 1672856277382934317688463)
    # btilde8 = convert(T, -74791376208282344108625901 // 737588957781464692067010)
    # btilde9 = convert(T, 1 // 30)

    # extra = Vern6ExtraStages(T, T2)
    # interp = Vern6InterpolationCoefficients(T)

    Vern6Tableau(
        c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63,
        a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, a91, a94,
        a95, a96, a97, a98,
        # btilde1, btilde4, btilde5, btilde6, btilde7, btilde8, btilde9, extra, interp
    )
end
