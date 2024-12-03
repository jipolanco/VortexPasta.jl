export ConstantFourierNormalFluidForcing,
       FourierNormalFluidForcing

"""
    FourierNormalFluidForcing <: NormalFluidForcing

Normal fluid forcing implemented in Fourier space.

See also [`NormalFluidForcing`](@ref).
"""
abstract type FourierNormalFluidForcing <: AbstractForcing end
