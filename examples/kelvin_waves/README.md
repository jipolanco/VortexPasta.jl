# Simulation of unforced and forced Kelvin waves

## Running the example

From this directory, open Julia with:

```bash
julia [-t NTHREADS] --project=.
```

where the `-t NTHREADS` is optional and allows using multithreading.

Then, install required packages by entering package mode (with `]`):

```julia
(kelvin_waves) pkg> instantiate
```

In particular, this will install the latest release of VortexPasta.
If you want to use the version in this repository instead, do:

```julia
(kelvin_waves) pkg> dev ../..
```

Finally, the example scripts can be executed using `include` (or, even better,
using `includet` from the [Revise.jl](https://github.com/timholy/Revise.jl) package):

```julia
includet("forced_kelvin_waves.jl")
```

If using Revise, one can make modifications to the `run_forced_lines` or other
called functions, and these modifications will be automatically taken into
account in next calls to `run_forced_lines` from the same Julia session.
