# Parallelisation and clusters

VortexPasta.jl can take advantage of thread-based CPU parallelisation (equivalent to
OpenMP).
In particular, it can be run on a single node of a computing cluster.

## Starting Julia with multiple threads

There are [a few ways](https://docs.julialang.org/en/v1/manual/multi-threading/) of starting
Julia with multiple threads:

1. either via the [`JULIA_NUM_THREADS`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_NUM_THREADS) environment variable,

2. or via the `-t / --threads` command-line option.

If both are used, the second takes precedence.

The command-line flag can be used in two ways:

```bash
$ julia -t 8     # start Julia with 8 threads
$ julia -t auto  # start Julia with the number of CPUs available to this Julia process
```

The second option can be useful in particular when running [SLURM jobs](@ref
slurm-job), as Julia will use the number of CPUs associated to the SLURM allocation.


## Pinning threads

When Julia is started with multiple threads, it can (and often does) assign
more than one thread to the same CPU, even when the requested number of threads
is smaller or equal to the total number of CPUs.
This is clearly suboptimal.

The [ThreadPinning.jl](https://github.com/carstenbauer/ThreadPinning.jl) package
solves this issue.
The easiest way to use it by putting [the following lines](https://carstenbauer.github.io/ThreadPinning.jl/stable/examples/ex_pinning_julia_threads/#pinthreads) at the top of your Julia script:

```julia
using ThreadPinning
pinthreads(:cores)
```

Note that this should be changed [when using SLURM](@ref Pinning-SLURM-threads).

One can then check that threads are correctly pinned to separate CPUs using [`threadinfo`](https://carstenbauer.github.io/ThreadPinning.jl/stable/examples/ex_pinning_julia_threads/#threadinfo_example).

## [Running SLURM jobs](@id slurm-job)

### [Submitting jobs](@id slurm-submitting-jobs)

Here is a sample SLURM script for running a simulation on
which can be submitted using `sbatch`:

```bash
#!/bin/bash

#SBATCH --job-name="JOB_NAME"
#SBATCH --partition=PARTITION_NAME_ON_CLUSTER
#SBATCH --time=1:00:00
#SBATCH --distribution=block:block
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --hint=nomultithread
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --mem=120G

echo " - SLURM_JOB_ID = $SLURM_JOB_ID"
echo " - SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo " - SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo " - SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"

srun --cpu-bind=verbose,cores \
  --cpus-per-task=$SLURM_CPUS_PER_TASK \
  --ntasks-per-node=$SLURM_NTASKS_PER_NODE \
  --threads-per-core=$SLURM_THREADS_PER_CORE \
  julia -t auto --heap-size-hint=100G --project=. script.jl
```

Here the number of threads is set via the `--cpus-per-task` option.

Some other notes:

- the precise `SBATCH` flags that need to be passed depend on the actual cluster.
  For example, in some clusters one may also need to pass `--constraint` and `--account` flags,
  while the `--partition` flag is not necessarily needed.

- the `--exclusive` flag is optional; it is recommended if one wants to use
  a full node (i.e. if the number of requested CPUs corresponds to the total
  number of CPUs on a node);

- above, we passed the `--heap-size-hint=100G` option to the `julia` command.
  This may help avoid out of memory errors, by telling Julia's garbage
  collector that the total used memory should stay below the specified value.
  Just to be sure, we also explicitly requested a (slightly larger) amount
  of memory to SLURM using the `--mem` option;

- the `--hint=nomultithread` option tells SLURM to avoid using hyperthreading,
  which is generally not good for performance.

### Pinning SLURM threads

As mentioned [above](@ref Pinning-threads), it is a good idea to pin Julia threads
to the CPUs available to the Julia process.

When using SLURM, one can achieve this by using the `:affinitymask` criterion in ThreadPinning's [`pinthreads`](https://carstenbauer.github.io/ThreadPinning.jl/stable/refs/api_pinning/#ThreadPinning.pinthreads):

```julia
using ThreadPinning
pinthreads(:affinitymask)
```

It can be convenient to have a Julia script which does the Right Thing (TM)
depending on whether it runs within a SLURM job or not.
To achieve this, one can do:

```julia
using ThreadPinning

if haskey(ENV, "SLURM_JOB_ID")
    pinthreads(:affinitymask)
else
    pinthreads(:cores)
end
```
## Using MKL FFT routines

The default [`NonuniformFFTsBackend`](@ref) in VortexPasta.jl computes threaded
FFTs using the [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) package, which
by default wraps the FFTW libraries written in C.
However, FFTW.jl has an [unresolved
issue](https://github.com/JuliaMath/FFTW.jl/issues/236) which can be
encountered (somewhat randomly) when computing FFTs using a large number of
threads.

### Switching to MKL in FFTW.jl

One workaround is to switch to the FFT implementation in Intel's MKL libraries, which don't seem to display this issue.
The FFTW.jl package makes it easy to [switch](https://github.com/JuliaMath/FFTW.jl?tab=readme-ov-file#mkl) to the MKL implementation via their [FFTW interface](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/fftw3-interface-to-onemkl.html).

One simply needs to do:

```julia
using FFTW
FFTW.set_provider!("mkl")
```

and restart Julia. This will create (or update) a `LocalPreferences.toml` file next to the `Project.toml` file associated to the active Julia project.

### Correctly using threads with MKL

The above change is not enough if one wants MKL's FFTs to be efficient when
using threads.
One also needs to set the following environment variables (for example in a [SLURM script](@ref slurm-submitting-jobs)):

```bash
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # on SLURM
export MKL_NUM_THREADS=$NUMBER_OF_JULIA_THREADS  # in general
export MKL_DYNAMIC=false
```

The `MKL_DYNAMIC=false` option tells MKL [not to mess around with thread pinning](https://carstenbauer.github.io/ThreadPinning.jl/stable/examples/ex_blas/#Intel-MKL).

Secondly, one also needs to add:

```julia
using MKL
```

at the start of the Julia script to be run (one may need to `]add MKL` first).
Failing to do this can really degrade performance.

Note that one can also set the environment variables directly in the Julia script.
Including [thread pinning](@ref Pinning-SLURM-threads), the beginning of the
Julia script could look like:

```julia
using MKL
using ThreadPinning

ENV["MKL_NUM_THREADS"] = Threads.nthreads()  # same as number of Julia threads
ENV["MKL_DYNAMIC"] = false

if haskey(ENV, "SLURM_JOB_ID")
    pinthreads(:affinitymask)
else
    pinthreads(:cores)
end

# Tell threadinfo to give us information about BLAS (and MKL) and optionally about the SLURM set-up.
threadinfo(blas = true, hints = true, slurm = haskey(ENV, "SLURM_JOB_ID"))
```

Note that, with the `hints = true` option, ThreadPinning will complain about
our choice of using `MKL_NUM_THREADS = number_of_julia_threads`. This warning
can be ignored, since FFTs are executed from a single Julia thread and it's
therefore what we want.
