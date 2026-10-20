module Devices

using KernelAbstractions: KernelAbstractions as KA
export CPU, GPU, DeviceType

# Tools for selecting between CPU and GPU devices (via KernelAbstractions.jl)

struct CPU end
const GPU = KA.Backend
const DeviceType = Union{CPU, GPU}

to_ka_backend(backend::GPU) = backend
to_ka_backend(::CPU) = KA.CPU()  # this is POCLBackend() on KA ≥ 0.10

# Make our CPU backend compatible with KernelAbstractions
KA.device!(::CPU, device) = nothing  # not relevant on CPU
KA.device(::CPU) = 1
KA.ndevices(::CPU) = 1
KA.synchronize(::CPU) = nothing

# These use the KA.CPU backend, which is POCLBackend on KA ≥ 0.10.
KA.allocate(::CPU, ::Type{T}, args...) where {T} = KA.allocate(KA.CPU(), T, args...)::Array{T}
KA.zeros(::CPU, ::Type{T}, args...) where {T} = KA.zeros(KA.CPU(), T, args...)::Array{T}

end
