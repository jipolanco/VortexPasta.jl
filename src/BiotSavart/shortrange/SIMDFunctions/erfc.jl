# Ported from xsimd:
# https://github.com/xtensor-stack/xsimd/blob/4d185a6de75ff51fd466064c8d91557ae459105a/include/xsimd/arch/common/xsimd_common_math.hpp#L528
#
# /***************************************************************************
#  * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
#  * Martin Renou                                                             *
#  * Copyright (c) QuantStack                                                 *
#  * Copyright (c) Serge Guelton                                              *
#  *                                                                          *
#  * Distributed under the terms of the BSD 3-Clause License.                 *
#  *                                                                          *
#  * The full license is in the file LICENSE, distributed with this software. *
#  ****************************************************************************/

# computes erf(a0)/a0
# x is sqr(a0) and 0 <= abs(a0) <= 2/3
@inline function erf1(::Type{Float32}, x)
    a = (
        reinterpret(Float32, 0x3f906eba), #   1.128379154774254e+00
        reinterpret(Float32, 0xbec0937e), #  -3.761252839094832e-01
        reinterpret(Float32, 0x3de70f22), #   1.128218315189123e-01
        reinterpret(Float32, 0xbcdb61f4), #  -2.678010670585737e-02
        reinterpret(Float32, 0x3ba4468d), #   5.013293006147870e-03
        reinterpret(Float32, 0xba1fc83b), #  -6.095205117313012e-04
    )
    evalpoly(x, a)
end

# computes erfc(x)*exp(sqr(x))
# x >=  2/3
@inline function erfc2(::Type{Float32}, x)
    a = (
        reinterpret(Float32, 0x3f0a0e8b), #   5.392844046572836e-01
        reinterpret(Float32, 0xbf918a62), #  -1.137035586823118e+00
        reinterpret(Float32, 0x3e243828), #   1.603704761054187e-01
        reinterpret(Float32, 0x3ec4ca6e), #   3.843569094305250e-01
        reinterpret(Float32, 0x3e1175c7), #   1.420508523645926e-01
        reinterpret(Float32, 0x3e2006f0), #   1.562764709849380e-01
        reinterpret(Float32, 0xbfaea865), #  -1.364514006347145e+00
        reinterpret(Float32, 0x4050b063), #   3.260765682222576e+00
        reinterpret(Float32, 0xc0cd1a85), #  -6.409487379234005e+00
        reinterpret(Float32, 0x40d67e3b), #   6.702908785399893e+00
        reinterpret(Float32, 0xc0283611), #  -2.628299919293280e+00
    )
    evalpoly(x, a)
end

@inline function erfc3(::Type{Float32}, x)
    a = (
        reinterpret(Float32, 0x3f7ffffe), #   9.9999988e-01
        reinterpret(Float32, 0xbe036d7e), #  -1.2834737e-01
        reinterpret(Float32, 0xbfa11698), #  -1.2585020e+00
        reinterpret(Float32, 0xbffc9284), #  -1.9732213e+00
        reinterpret(Float32, 0xc016c985), #  -2.3560498e+00
        reinterpret(Float32, 0x3f2cff3b), #   6.7576951e-01
        reinterpret(Float32, 0xc010d956), #  -2.2632651e+00
        reinterpret(Float32, 0x401b5680), #   2.4271545e+00
        reinterpret(Float32, 0x41aa8e55), #   2.1319498e+01
    )
    (one(x) - x) * evalpoly(x, a)
end

# computes erf(a0)/a0
# x is sqr(a0) and 0 <= abs(a0) <= 0.65
@inline function erf1(::Type{Float64}, x)
    a = (
        reinterpret(Float64, 0x3ff20dd750429b61), # 1.12837916709551
        reinterpret(Float64, 0x3fc16500f106c0a5), # 0.135894887627278
        reinterpret(Float64, 0x3fa4a59a4f02579c), # 4.03259488531795E-02
        reinterpret(Float64, 0x3f53b7664358865a), # 1.20339380863079E-03
        reinterpret(Float64, 0x3f110512d5b20332), # 6.49254556481904E-05
    )
    b = (
        reinterpret(Float64, 0x3ff0000000000000), # 1
        reinterpret(Float64, 0x3fdd0a84eb1ca867), # 0.453767041780003
        reinterpret(Float64, 0x3fb64536ca92ea2f), # 8.69936222615386E-02
        reinterpret(Float64, 0x3f8166f75999dbd1), # 8.49717371168693E-03
        reinterpret(Float64, 0x3f37ea4332348252), # 3.64915280629351E-04
    )
    ax = evalpoly(x, a)
    bx = evalpoly(x, b)
    Base.FastMath.div_fast(ax, bx)
end

# computes erfc(x)*exp(x*x)
# 0.65 <= abs(x) <= 2.2
@inline function erfc2(::Type{Float64}, x)
    a = (
        reinterpret(Float64, 0x3feffffffbbb552b), # 0.999999992049799
        reinterpret(Float64, 0x3ff54dfe9b258a60), # 1.33154163936765
        reinterpret(Float64, 0x3fec1986509e687b), # 0.878115804155882
        reinterpret(Float64, 0x3fd53dd7a67c7e9f), # 0.331899559578213
        reinterpret(Float64, 0x3fb2488a6b5cb5e5), # 7.14193832506776E-02
        reinterpret(Float64, 0x3f7cf4cfe0aacbb4), # 7.06940843763253E-03
        reinterpret(Float64, 0x0000000000000000), # 0
    )
    b = (
        reinterpret(Float64, 0x3ff0000000000000), # 1
        reinterpret(Float64, 0x4003adeae79b9708), # 2.45992070144246
        reinterpret(Float64, 0x40053b1052dca8bd), # 2.65383972869776
        reinterpret(Float64, 0x3ff9e677c2777c3c), # 1.61876655543871
        reinterpret(Float64, 0x3fe307622fcff772), # 0.594651311286482
        reinterpret(Float64, 0x3fc033c113a7deee), # 0.126579413030178
        reinterpret(Float64, 0x3f89a996639b0d00), # 1.25304936549413E-02
    )
    ax = evalpoly(x, a)
    bx = evalpoly(x, b)
    Base.FastMath.div_fast(ax, bx)
end

# computes erfc(x)*exp(x*x)
# 2.2 <= abs(x) <= 6
@inline function erfc3(::Type{Float64}, x)
    a = (
        reinterpret(Float64, 0x3fefff5a9e697ae2), # 0.99992114009714
        reinterpret(Float64, 0x3ff9fa202deb88e5), # 1.62356584489367
        reinterpret(Float64, 0x3ff44744306832ae), # 1.26739901455873
        reinterpret(Float64, 0x3fe29be1cff90d94), # 0.581528574177741
        reinterpret(Float64, 0x3fc42210f88b9d43), # 0.157289620742839
        reinterpret(Float64, 0x3f971d0907ea7a92), # 2.25716982919218E-02
        reinterpret(Float64, 0x0000000000000000), # 0
    )
    b = (
        reinterpret(Float64, 0x3ff0000000000000), # 1
        reinterpret(Float64, 0x400602f24bf3fdb6), # 2.75143870676376
        reinterpret(Float64, 0x400afd487397568f), # 3.37367334657285
        reinterpret(Float64, 0x400315ffdfd5ce91), # 2.38574194785344
        reinterpret(Float64, 0x3ff0cfd4cb6cde9f), # 1.05074004614827
        reinterpret(Float64, 0x3fd1d7ab774bb837), # 0.278788439273629
        reinterpret(Float64, 0x3fa47bd61bbb3843), # 4.00072964526861E-02
    )
    ax = evalpoly(x, a)
    bx = evalpoly(x, b)
    Base.FastMath.div_fast(ax, bx)
end

# computes erfc(rx)*exp(rx*rx)
# x >=  6 rx = 1/x
# -- this one is not needed?
# @inline function erfc4(::Type{Float64}, x)
#     a = (
#         reinterpret(Float64, 0xbc7e4ad1ec7d0000), # -2.627435221016534e-17
#         reinterpret(Float64, 0x3fe20dd750429a16), # 5.641895835477182e-01
#         reinterpret(Float64, 0x3db60000e984b501), # 2.000889609806154e-11
#         reinterpret(Float64, 0xbfd20dd753ae5dfd), # -2.820947949598745e-01
#         reinterpret(Float64, 0x3e907e71e046a820), # 2.457786367990903e-07
#         reinterpret(Float64, 0x3fdb1494cac06d39), # 4.231311779019112e-01
#         reinterpret(Float64, 0x3f34a451701654f1), # 3.149699042180451e-04
#         reinterpret(Float64, 0xbff105e6b8ef1a63), # -1.063940737150596e+00
#         reinterpret(Float64, 0x3fb505a857e9ccc8), # 8.211757799454056e-02
#         reinterpret(Float64, 0x40074fbabc514212), # 2.913930388669777e+00
#         reinterpret(Float64, 0x4015ac7631f7ac4f), # 5.418419628850713e+00
#         reinterpret(Float64, 0xc0457e03041e9d8b), # -4.298446704382794e+01
#         reinterpret(Float64, 0x4055803d26c4ec4f), # 8.600373238783617e+01
#         reinterpret(Float64, 0xc0505fce04ec4ec5), # -6.549694941594051e+01
#     )
#     evalpoly(x, a)
# end

@inline function _erfc(::Type{T}, self) where {T <: Float32}
    x = abs(self)
    test0 = self < zero(self)
    r1 = zero(self)
    test1 = 3f0 * x < 2f0
    z = Base.FastMath.div_fast(x, one(x) + x)
    if any(test1)
        r1 = erfc3(T, z)
        if all(test1)
            return vifelse(test0, 2f0 - r1, r1)
        end
    end
    z -= 0.4f0
    r2 = exp(-x * x) * erfc2(T, z)
    r1 = vifelse(test1, r1, r2)
    # r1 = vifelse(x == Inf32, zero(r1), r1)
    vifelse(test0, 2f0 - r1, r1)
end

@inline function _erfc(::Type{T}, self) where {T <: Float64}
    x = abs(self)
    xx = x * x
    test0 = self < 0.0
    test1 = x < 0.65
    r1 = zero(x)
    if any(test1)
        r1 = one(x) - x * erf1(T, xx)
        if all(test1)
            return vifelse(test0, 2.0 - r1, r1)
        end
    end
    test2 = x < 2.2
    test3 = test2 & ~test1
    ex = exp(-xx)
    if any(test3)
        z = ex * erfc2(T, x)
        r1 = vifelse(test1, r1, z)
        if all(test1 | test3)
            return vifelse(test0, 2.0 - r1, r1)
        end
    end
    z = ex * erfc3(T, x)
    r1 = vifelse(test2, r1, z)
    # r1 = vifelse(x == Inf, zero(r1), r1)
    vifelse(test0, 2.0 - r1, r1)
end

@inline function erfc(x)
    T = ftype(x)
    _erfc(T, x)
end
