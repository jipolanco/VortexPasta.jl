# VortexPasta.jl

A vortex filament solver for the simulation of quantum vortex flows.

## Introduction

[Quantum vortices](https://en.wikipedia.org/wiki/Quantum_vortex) are one of the most important features of superfluids such as low-temperature liquid helium or Bose--Einstein condensates.
These are *topological defects* -- which in three dimensions take the form of spatial curves -- where the velocity of the superfluid is singular.
The main property of quantum vortices is that they have a *quantised circulation*.
That is, they induce a rotating velocity field about them whose intensity is quantised in terms of the quantum of circulation ``κ = h/m`` where ``h`` is [Planck's constant](https://en.wikipedia.org/wiki/Planck_constant) and ``m`` is the mass of one boson (e.g. the mass of a ⁴He atom).

In the case of helium-4, superfluidity takes place at temperatures below about 2.17 K,
and the thickness of quantum vortices is about 1 Å (``10^{-10}`` m).
Their atomic-scale thickness justifies treating quantum vortices as infinitesimal lines when describing the macroscopic flow induced by them.

### The vortex filament model

The **vortex filament model** (VFM) describes the motion of thin vortex lines in three-dimensional space.
Each vortex line induces a velocity field about it given by the [Biot--Savart law](https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law):
```math
\bm{v}(\bm{x}) =
\frac{Γ}{4π} ∮_{\mathcal{C}} \frac{(\bm{s} - \bm{x}) \times \mathrm{d}\bm{s}}{|\bm{s} - \bm{x}|^3}
```
where ``Γ`` is the vortex circulation (equal to ``κ`` for quantum vortices),
and ``\bm{s}`` is a point along the vortex line.
The above equation derives from the vorticity field:
```math
\bm{ω}(\bm{x}) ≡ \bm{\nabla} × \bm{v}(\bm{x})
= Γ ∮_{\mathcal{C}} δ(\bm{s} - \bm{x}) \, \mathrm{d}\bm{s}
```
where ``δ`` is Dirac delta function.
That is, the vorticity field is singular and localised at the locations of quantum vortices.

The Biot--Savart law describes in particular the motion induced by vortex filaments on themselves and on surrounding vortex lines.
The VFM thus describes the collective motion of a set of mutually-interacting vortex filaments which obey the Biot--Savart law.
Note that the Biot--Savart integral is singular when evaluated at a vortex location ``\bm{s}' ∈ \mathcal{C}``, and the integral must be desingularised by taking into account the finite thickness of the vortex core.
The VFM also accounts for *vortex reconnections*, which occur when two vortex segments are sufficiently close to each other and which change the topology of the vortex system.

### Numerical considerations

From a numerical standpoint, simulating large sets of quantum vortices using the VFM is expensive.
Indeed, assuming that vortex lines are discretised by a set of ``N`` points in space, the VFM requires computing the velocity of each vortex point ``\bm{s}_i`` (for ``i ∈ [1, N]``) induced by the locations of *all* vortex points.
Note that the influence of a vortex on the surrounding space decays relatively slowly with the distance ``r`` (the induced velocity typically decays as ``1/r``), and thus one needs to account for all pair interactions -- even over long distances.
A naive implementation, where all such pair interactions are computed one by one, leads to a computational complexity of ``\mathcal{O}(N^2)``, that is, the computational cost increases quickly with the number of discretisation points (or vortices).
For this reason, simulating fully turbulent flows is prohibitively costly with naive methods.

The VortexPasta.jl solver provides an efficient and accurate implementation of the VFM, specifically intended for vortex flows in **periodic domains**.
Its efficiency comes from the use of a splitting technique, derived from the
[Ewald summation](https://en.wikipedia.org/wiki/Ewald_summation) method, which splits the Biot--Savart integral into a *short-range* and a *long-range* parts.
This kind of methods is routinely used in the context of molecular dynamics simulations to compute electrostatic Coulomb interactions between point charges.
In the modified short-range integral, pair interactions decay exponentially with ``r``, which means that one can effectively ignore interactions beyond some cut-off distance ``r_{\text{cut}}``.
Meanwhile, the long-range integral is smooth (non-singular) at ``r = 0``, and thus can be indirectly computed using fast Fourier transforms (FFTs), leading to a complexity of ``\mathcal{O}(N \log N)``.
Conveniently, the long-range velocity field is nothing but the velocity induced by
a coarse-grained vorticity field, which enables a simple interpretation of the long-range component, as well as the use of tools and quantities commonly used for classical flows (such as energy spectra, …).

## Installation

### Installing Julia

The easiest way of installing (and updating) Julia is using the official [`juliaup`](https://github.com/JuliaLang/juliaup) installer.

On Linux and Mac, one can install `juliaup` by executing in a terminal

```bash
curl -fsSL https://install.julialang.org | sh
```

This will automatically install the latest Julia release and make it available from your terminal using the `julia` command.

To get started with using Julia, the [Modern Julia Workflows blog](https://modernjuliaworkflows.github.io/) is a great source of practical information.
See in particular the section on [writing Julia code](https://modernjuliaworkflows.github.io/pages/writing/writing/).

In particular, it is highly recommended to install the [Revise.jl](https://github.com/timholy/Revise.jl/) package in the default environment.
This package keeps track of modifications of included code when working from a Julia session, which means that one doesn't need to re-include a Julia file each time it is modified.
To install it, launch Julia in a terminal, enter package mode using `]`, and then:

```julia-repl
(@v1.10) pkg> add Revise
```

Another very widely used package is [OhMyREPL.jl](https://github.com/KristofferC/OhMyREPL.jl), which adds syntax highlighting to the REPL among other nice things.
It can be installed in the same way as Revise.jl.

One usually wants to load these packages every time a Julia session is started.
This can be done automatically by creating a `~/.julia/config/startup.jl` with the following content:

```julia
try
    using Revise
    using OhMyREPL
catch e
    @warn "Error importing Revise or OhMyREPL. Use `]add Revise OhMyREPL` to install them."
end
```

### Local environments

Above, both packages have been installed in the default environment (`@v1.X`), which is convenient because we want to load these packages no matter which Julia project we're currently working on.

For more specific things (such as doing simulations using VortexPasta.jl), it is [recommended](https://modernjuliaworkflows.github.io/pages/writing/writing/#environments) to install things in a *local* environment, which allows (among other things) to avoid conflicts between package versions when working on different Julia projects.

To create a new local environment in the current folder, launch Julia using:[^1]

```bash
julia --project=.
```

Packages will then be installed in this environment.
For example, enter package mode using `]` and install the [GLMakie.jl](https://docs.makie.org/stable/) package:

```julia
(simulation) pkg> add GLMakie
```

This will install GLMakie package (for plotting) and all its dependencies.
It will also precompile installed packages, which can take a few minutes (especially for GLMakie, which is a large package with many dependencies), but which makes loading the package faster when it is actually imported from a Julia session or script.

One can check that the project now contains a single package:

```julia
(simulation) pkg> status
Status `/path/to/simulation/Project.toml`
  [e9467ef8] GLMakie v0.8.9
```

Moreover, the current directory now contains two new files, `Project.toml` and `Manifest.toml`, which describe the packages installed by local environment.
The `Project.toml` file lists the packages that we have directly installed (so just GLMakie for now).
The `Manifest.toml` file is much larger, as it also contains the list of *all* installed packages associated to the environment, including the direct and indirect dependencies of GLMakie, and the precise versions of the installed packages.
This file is automatically generated and is very useful for reproducibility.

Note that, once the new project has been created, one can simply start Julia using

```bash
julia --project
```

from the same directory (`simulation` in this example).
Julia will detect the presence of the `Project.toml` file and automatically enable the local environment.

### Installing VortexPasta.jl

It is recommended to install VortexPasta.jl in a local environment, for example in the same environment as in the previous example (which already includes GLMakie).
Since it is not currently registered in a Julia package registry, it should be installed as:

```julia-repl
(simulation) pkg> add git@gitlab.in2p3.fr:jipolanco/VortexPasta.jl.git
```

As before, this will install all dependencies and precompile VortexPasta.jl and other installed packages.
Now the status of the local environment should look something like:

```julia-repl
(simulation) pkg> status
Status `~/Work/Projects/Vortices/Codes/VortexPasta/docs/simulation/Project.toml`
  [e9467ef8] GLMakie v0.8.9
  [3d0f1e53] VortexPasta v0.1.0 `git@gitlab.in2p3.fr:jipolanco/VortexPasta.jl.git#master`
```

!!! warning

    If the `add` command asks for an SSH key location and then fails with a `GitError`, [you may need to set](https://github.com/GunnarFarneback/LocalRegistry.jl/blob/master/docs/ssh_keys.md#2-using-an-external-git-binary-with-julias-package-manager) the [`JULIA_PKG_USE_CLI_GIT`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_PKG_USE_CLI_GIT) environment variable to `true`.
    For instance, add the line
    ```julia
    ENV["JULIA_PKG_USE_CLI_GIT"] = "true"
    ```
    to your `~/.julia/config/startup.jl` file and then restart Julia.


## Running simulations and analysis

TODO: refer to the first tutorial?


[^1]: Alternatively, one can start Julia *without* the `--project` flag, and then `] activate .` in package mode, as described in the [official docs](https://pkgdocs.julialang.org/dev/getting-started/#Getting-Started-with-Environments).
