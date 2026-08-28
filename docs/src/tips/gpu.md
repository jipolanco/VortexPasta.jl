# GPU usage

Currently, it is possible to accelerate the computation of both [short-range](@ref Short-range-velocity) and [long-range
interactions](@ref Long-range-velocity) using GPUs.
The respective GPU-compatible backends are [`CellListsBackend`](@ref) and [`NonuniformFFTsBackend`](@ref),
see their documentations for details.

## Using a local CUDA toolkit

By default, CUDA.jl will download the CUDA toolkit when it is first used.
One may want to avoid it and instead use a CUDA toolkit already installed in the system.
This is particularly important on **HPC clusters**, where an internet connection may not be available on compute nodes (GPU-enabled) but only on login nodes (without GPUs).

To do this, one should first determine the version of the available CUDA toolkit.
This may be obtained by running `nvcc --version` on the system where computations will be performed (e.g. from a compute node).
Look for a line similar to `Cuda compilation tools, release 13.2, V13.2.78` (here the version is `13.2`).
On HPC clusters, one may first need to load a CUDA module, e.g. `module load cuda` (but this will depend on the cluster).

Then, the configuration should be performed in the two steps detailed below.
The first one, done in the global environment, should only be done once (until
either the CUDA or the Julia version change, in which case it should be repeated).
The second step is simpler and should be done each time CUDA.jl is to be used in a local project (e.g. along VortexPasta.jl).

### 1. Installing CUDA.jl on the global environment

Now, install CUDA.jl on the _global_ Julia environment[^1] from a machine with internet access (typically the login node, without GPU).
For this, launch `julia` _without_ the `--project` flag (to use the global environment), and
then:

```julia-repl
julia> using Pkg; Pkg.add("CUDA")  # install CUDA.jl

julia> using CUDA

julia> CUDA.set_runtime_version!(v"13.2"; local_toolkit = true)
```

The `13.2` should be replaced with the CUDA version found using `nvcc --version`.
See the [CUDA.jl docs](https://cuda.juliagpu.org/dev/installation/overview/#Using-a-local-CUDA) for more details.

This should generate a `LocalPreferences.toml` file under `$JULIA_DEPOT_PATH/environments/v1.12/` (replace `v1.12` with the current Julia version), which should look as follows:

```toml
[CUDA_Compiler_jll]
local = "true"
version = "13.2"

[CUDA_Runtime_jll]
local = "true"
version = "13.2"
```

Afterwards, it may be helpful to launch `julia` again and run:

```julia-repl
julia> using CUDA

julia> CUDA.precompile_runtime()
```

### 2. Using CUDA.jl in a local environment

After CUDA.jl has been configured in the global environment, one would want to
reuse the same configuration in a local environment (which may contain other
packages such as VortexPasta.jl to be used in a local project).

To do this, launch Julia from the environment associated to the local project
by doing something like `julia --project=.` (see [Local environments](@ref
julia-local-environments) for more details).

Then, "install" CUDA.jl in this environment:

```julia-repl
julia> using Pkg; Pkg.add("CUDA")
```

That's it!

This local installation seems to be needed even if CUDA.jl is already in the
global environment, probably because the versions of other packages in the
local project are compatible with the installed CUDA.jl version (so they might
be downgraded when adding CUDA).

Finally, one can check the CUDA.jl configuration on a GPU-enabled node:

```julia-repl
julia> using CUDA

julia> CUDA.versioninfo()
CUDA toolchain: 
- runtime 13.2.0, local installation
- driver 595.71.5 for 13.3
- compiler 13.2.78, local installation

CUDA libraries: 
- cuBLAS: 13.4.0
- cuSPARSE: 12.7.10
- cuSOLVER: 12.2.0
- cuFFT: 12.2.0
- cuRAND: 10.4.2
- CUPTI: 2026.1.1 (API 13.2.1)
- NVML: 13.0.0+595.71.5

Julia packages: 
- CUDACore: 6.3.0
- GPUArrays: 11.5.13
- GPUCompiler: 2.4.2
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.1+0
- CUDA_Compiler_jll: 0.5.1+0
- CUDA_Runtime_jll: 0.24.1+0
- CUDA_Runtime_Discovery: 2.1.0
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.7
- LLVM: 18.1.7

Preferences:
- CUDA_Runtime_jll.version: 13.2
- CUDA_Runtime_jll.local: true
- CUDA_Compiler_jll.version: 13.2
- CUDA_Compiler_jll.local: true

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 79.177 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.2 (LLVM: sm_90a / PTX 9.0)
```

One can verify that the `Preferences` section contains the same values as the `LocalPreferences.toml` file generated in the global environment.

[^1]: It's not strictly necessary to do this step in the global environment (it could be done each time CUDA.jl is installed on a local environment), but it's a convenient way of doing this once and for all.
