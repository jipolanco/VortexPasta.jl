# Vectorised implementation of erf(x) adapted from VectorizationBase.jl
# See https://github.com/JuliaSIMD/VectorizationBase.jl/blob/129a0e533202ee3bf1b60a925f74ad153e8bcdd7/src/special/verf.jl.

module SIMDFunctions

using SIMD: Vec, vifelse

@inline vfmadd(a, b, c) = muladd(a, b, c)
@inline vfnmadd(a, b, c) = vfmadd(-a, b, c)
@inline vfmsub(a, b, c) = vfmadd(a, b, -c)
@inline vfnmsub(a, b, c) = -vfmadd(a, b, c)

@inline function erf_kernel_l9(x::Union{Float64, Vec{<:Any, Float64}})
    y = vfmadd(x, 6.49254556481904e-5, 0.00120339380863079)
    z = vfmadd(x, 0.000364915280629351, 0.00849717371168693)
    y = vfmadd(x, y, 0.0403259488531795)
    z = vfmadd(x, z, 0.0869936222615386)
    y = vfmadd(x, y, 0.135894887627278)
    z = vfmadd(x, z, 0.453767041780003)
    y = vfmadd(x, y, 1.12837916709551)
    z = vfmadd(x, z, 1.0)
    res = y / z
    res + vfnmadd(
        8.13035732583580548119176214095147680242299108680658322550212199643608432312984e-16,
        x,
        2.679870718713541577762315771546743935213688171443421284936882991116060314650674e-15
    )
    # Base.FastMath.div_fast(w * y, z)
end

# vscalef(b::Bool, x, y, z) = b ? (x * (2^y)) : z
vscalef(b::Vec, x, y, z) = vifelse(b, x * (2^y), z)

@inline function erf_v56f_kernel(v3f, v4f)
    v29f = v4f * -1.4426950408889634
    v30f = round(v29f)
    v31f = vfmsub(v30f, -0.6931471803691238, v4f)
    v32f = v30f * 1.9082149292705877e-10
    v33f = v31f - v32f
    v34f = v33f * v33f
    v35f = vfmadd(v34f, 4.1381367970572385e-8, -1.6533902205465252e-6)
    v36f = vfmadd(v34f, v35f, 6.613756321437934e-5)
    v37f = vfmadd(v34f, v36f, -0.0027777777777015593)
    v39f = vfnmsub(v34f, v37f, 0.16666666666666602)
    v40f = vfmadd(v34f, v39f, v33f)
    v41f = v40f * v33f
    v42f = 2.0 - v40f
    v43f = Base.FastMath.div_fast(v41f, v42f)
    v44f = 1.0 - v32f
    v45f = v44f + v31f
    v46f = v45f + v43f
    m47 = v4f ≤ 708.3964185322641
    v54f = vscalef(m47, v46f, v30f, zero(v3f))
    # v56f = ifelse(v4f < 709.782712893384, v54f, Inf) # TODO: NaN should return v54f
end

@inline function erf_v97f_kernel(v3f)
    v91f = vfmadd(v3f, 0.0400072964526861, 0.278788439273629)
    v86f = vfmadd(v3f, 0.0225716982919218, 0.157289620742839)
    v92f = vfmadd(v3f, v91f, 1.05074004614827)
    v87f = vfmadd(v3f, v86f, 0.581528574177741)
    v93f = vfmadd(v3f, v92f, 2.38574194785344)
    v88f = vfmadd(v3f, v87f, 1.26739901455873)
    v94f = vfmadd(v3f, v93f, 3.37367334657285)
    v89f = vfmadd(v3f, v88f, 1.62356584489367)
    v95f = vfmadd(v3f, v94f, 2.75143870676376)
    v90f = vfmadd(v3f, v89f, 0.99992114009714)
    v96f = vfmadd(v3f, v95f, 1.0)
    v97f = v90f / v96f
end

@inline function erf_v71_kernel(v3f)
    v65f = vfmadd(v3f, 0.0125304936549413, 0.126579413030178)
    v60f = vfmadd(v3f, 0.00706940843763253, 0.0714193832506776)
    v66f = vfmadd(v3f, v65f, 0.594651311286482)
    v61f = vfmadd(v3f, v60f, 0.331899559578213)
    v67f = vfmadd(v3f, v66f, 1.61876655543871)
    v62f = vfmadd(v3f, v61f, 0.878115804155882)
    v68f = vfmadd(v3f, v67f, 2.65383972869776)
    v63f = vfmadd(v3f, v62f, 1.33154163936765)
    v69f = vfmadd(v3f, v68f, 2.45992070144246)
    v64f = vfmadd(v3f, v63f, 0.999999992049799)
    v70f = vfmadd(v3f, v69f, 1.0)
    v64f, v70f
end

@inline function verf(v0f::Vec{<:Any, Float64})
    v3f = abs(v0f)
    v4f = v0f * v0f
    m6 = v3f < 0.675
    if any(m6)
        v19f = v0f * erf_kernel_l9(v4f)
        all(m6) && return v19f
    else
        v19f = zero(v0f)
        # v19f = _vundef(v0f)
    end
    m23 = v3f < 2.725
    v56f = erf_v56f_kernel(v3f, v4f)
    if any(m23 & (~m6)) # any(0.675 < v3f < 2.2) # l58
        v64f, v70f = erf_v71_kernel(v3f)
        v71f = Base.FastMath.div_fast(v56f, v70f)
        v73f = vfnmadd(v71f, v64f, 1.0)
        v78f = copysign(v73f, v0f)
        v84f = vifelse(m6, v19f, v78f)
        all(m23) && return v84f
    else
        v84f = v19f
    end
    v97f = erf_v97f_kernel(v3f)
    v99f = vfnmadd(v97f, v56f, 1.0)
    # v99f = sub_ieee(1.0, v98f)
    v104f = copysign(v99f, v0f)
    v105f = vifelse(m23, v84f, v104f)
    v105f
end

@inline function verf(v0f::Vec{<:Any, Float32})
    v3f = abs(v0f)
    m4 = v3f < 0.6666667f0
    v8f = v3f * v3f
    if any(m4) # L7
        v9f = vfmadd(v8f, -0.00060952053f0, 0.005013293f0)
        v10f = vfmadd(v8f, v9f, -0.026780106f0)
        v11f = vfmadd(v8f, v10f, 0.11282183f0)
        v12f = vfmadd(v8f, v11f, -0.37612528f0)
        v13f = vfmadd(v8f, v12f, 1.1283791f0)
        v17f = v13f * v0f
        all(m4) && return v17f
    else
        v17f = zero(v0f)
        # v17f = _vundef(v0f)
    end
    v18f = v3f + 1.0f0
    v19f = Base.FastMath.div_fast(v3f, v18f)
    v20f = v19f - 0.4f0
    v23f = -1.442695f0 * v8f
    v24f = round(v23f)
    v25f = vfmsub(v24f, -0.6933594f0, v8f)
    v26f = vfmadd(v24f, 0.00021219444f0, v25f)
    v27f = vfmadd(v26f, 0.0013997796f0, 0.008395563f0)
    v28f = vfmadd(v26f, v27f, 0.0416654f0)
    v31f = v26f * v26f
    v29f = vfmadd(v26f, v28f, 0.16666277f0)
    v30f = vfmadd(v26f, v29f, 0.5f0)
    v32f = vfmadd(v30f, v31f, v26f)
    v33f = v32f + 1.0f0
    m34 = v8f ≤ 88.37626f0
    v42f = vscalef(m34, v33f, v24f, zero(v0f))
    v43f = vfmadd(v20f, -2.6283f0, 6.702909f0)
    v44f = vfmadd(v20f, v43f, -6.4094872f0)
    v45f = vfmadd(v20f, v44f, 3.2607658f0)
    v46f = vfmadd(v20f, v45f, -1.364514f0)
    v47f = vfmadd(v20f, v46f, 0.15627646f0)
    v48f = vfmadd(v20f, v47f, 0.14205085f0)
    v49f = vfmadd(v20f, v48f, 0.38435692f0)
    v50f = vfmadd(v20f, v49f, 0.16037047f0)
    v51f = vfmadd(v20f, v50f, -1.1370356f0)
    v52f = vfmadd(v20f, v51f, 0.5392844f0)
    v53f = v52f * v42f
    v54f = 1.0f0 - v53f
    m55 = v0f < zero(v0f)
    v56f = -v54f
    v57f = vifelse(m55, v56f, v54f)
    v58f = vifelse(m4, v17f, v57f)
    m59 = v3f ≠ Inf32
    m60 = v0f > zero(v0f)
    v61f = vifelse(m60, one(v0f), zero(v0f))
    v62f = vifelse(m55, one(v0f), zero(v0f))
    v63f = v61f - v62f
    v65f = vifelse(v0f == v0f, v63f, NaN32)
    v66f = vifelse(m59, v58f, v65f)
    v66f
end

end
