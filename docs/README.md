# Documentation

## Generating the documentation

To generate the documentation locally in your system, launch Julia from this folder with the `--project` flag:

```bash
julia --project
```

Then tap `]` to enter package mode and run the following commands:

```julia-repl
(docs) pkg> dev ..
(docs) pkg> instantiate
```

The first line will ensure that we're generating the documentation from the local version of this package (in the parent directory `..`).
The second command will make sure that all dependencies for generating the documentation are installed.

Finally, exit package mode (using backspace) and run the [`make.jl`](docs/make.jl) script:

```julia-repl
julia> include("make.jl")
```

If everything goes well, you can now open the `build/index.html` file in a browser to see the generated docs.
