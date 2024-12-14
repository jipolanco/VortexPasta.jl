# GPU usage

Currently, it is possible to accelerate the computation of [long-range
interactions](@ref Long-range-velocity) using GPUs.
To do this, one should set `backend_long = NonuniformFFTsBackend(...)`
when creating a [`ParamsBiotSavart`](@ref) and use the wanted GPU device as detailed in
[`NonuniformFFTsBackend`](@ref).

## Using a local CUDA toolkit

By default, CUDA.jl will download the CUDA toolkit when it is first used.
One may want to avoid it and instead use a CUDA toolkit already installed in the system.
This is particularly important on **HPC clusters**, where an internet connection may not be available on compute nodes (GPU-enabled) but only on login nodes (without GPUs).

To do this, one should first determine the version of the available CUDA toolkit.
This may be obtained by running `nvidia-smi` on the system where computations will be performed (e.g. from a compute node).
Look for something like `CUDA Version: 12.4` on the top.
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

julia> CUDA.set_runtime_version!(v"12.4"; local_toolkit = true)
```

The `12.4` should be replaced with the CUDA version found using `nvidia-smi`.
See the [CUDA.jl docs](https://cuda.juliagpu.org/dev/installation/overview/#Using-a-local-CUDA) for more details.

This should generate a `LocalPreferences.toml` file under `$JULIA_DEPOT_PATH/environments/v1.11/` (replace `v1.11` with the current Julia version), which should look as follows:

```toml
[CUDA_Runtime_jll]
local = "true"
version = "12.4"
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
CUDA runtime 12.4, local installation
CUDA driver 12.6
NVIDIA driver 550.90.7

CUDA libraries:
- CUBLAS: 12.4.5
- CURAND: 10.3.5
- CUFFT: 11.2.1
- CUSOLVER: 11.6.1
- CUSPARSE: 12.3.1
- CUPTI: 2024.1.1 (API 22.0.0)
- NVML: 12.0.0+550.90.7

Julia packages:
- CUDA: 5.5.2
- CUDA_Driver_jll: 0.10.4+0
- CUDA_Runtime_jll: 0.15.5+0
- CUDA_Runtime_Discovery: 0.3.5

Toolchain:
- Julia: 1.11.2
- LLVM: 16.0.6

Preferences:
- CUDA_Runtime_jll.version: 12.4
- CUDA_Runtime_jll.local: true

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 79.092 GiB / 79.647 GiB available)
```

One can verify that the `Preferences` section contains the same values as the `LocalPreferences.toml` file generated in the global environment.

[^1]: It's not strictly necessary to do this step in the global environment (it could be done each time CUDA.jl is installed on a local environment), but it's a convenient way of doing this once and for all.
