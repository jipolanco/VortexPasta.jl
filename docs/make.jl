using Documenter
using Documenter: Remotes
using DocumenterCitations
using Literate: Literate

using VortexPasta

# We load a few submodules to make @ref's work in Literate-generated files.
using VortexPasta.Filaments
using VortexPasta.BiotSavart

# For some reason we need to explicitly load these packages here to avoid failures with
# doctests and @example blocks on Gitlab.
using Rotations
using SpecialFunctions
using GLMakie
using CairoMakie

DocMeta.setdocmeta!(
    VortexPasta.PaddedArrays,
    :DocTestSetup,
    :(using VortexPasta.PaddedArrays),
)

DocMeta.setdocmeta!(
    VortexPasta.PredefinedCurves,
    :DocTestSetup,
    quote
        using StaticArrays: SVector
        using Rotations
        using VortexPasta.PredefinedCurves
    end,
)

DocMeta.setdocmeta!(
    VortexPasta.BasicTypes,
    :DocTestSetup,
    :(using VortexPasta.BasicTypes),
)

DocMeta.setdocmeta!(
    VortexPasta.Filaments,
    :DocTestSetup,
    quote
        using StaticArrays: SVector
        using VortexPasta: VortexPasta
        using VortexPasta.Filaments
    end,
)

# doctest(VortexPasta; fix = true)

struct Gitlab <: Remotes.Remote
    url :: String
end

Remotes.repourl(remote::Gitlab) = remote.url

# Example:
# https://gitlab.in2p3.fr/jipolanco/VortexPasta.jl/-/blob/master/src/Filaments/integrate.jl#L23-35
function Remotes.fileurl(remote::Gitlab, ref, filename, linerange)
    io = IOBuffer()
    print(io, Remotes.repourl(remote), "/-/blob/", ref, '/', filename)
    if linerange !== nothing
        a, b = first(linerange), last(linerange)
        print(io, "#L", a)
        if a != b
            print(io, "-", b)
        end
    end
    String(take!(io))
end

repo = Gitlab("https://gitlab.in2p3.fr/jipolanco/VortexPasta.jl")

##

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "biblio.bib");
    style = :authoryear,
)

Literate.markdown(
    joinpath(@__DIR__, "literate", "tutorials", "01-vortex_ring.jl"),
    joinpath(@__DIR__, "src", "tutorials");
    documenter = true,
)

makedocs(;
    sitename = "VortexPasta",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = Remotes.repourl(repo),
        edit_link = "master",
        size_threshold_ignore = [
            # "tutorials/01-vortex_ring.md",
            "modules/Filaments.md",  # this page has too much content so it's relatively large
        ],
        mathengine = KaTeX(),
    ),
    modules = [VortexPasta],
    pages = [
        "index.md",
        "Numerical considerations" => [
            "methods/Ewald.md",
        ],
        "Tutorials" => [
            "tutorials/01-vortex_ring.md",
        ],
        "Modules" => [
            "modules/PaddedArrays.md",
            "modules/PredefinedCurves.md",
            "modules/CellLists.md",
            "modules/BasicTypes.md",
            "modules/Quadratures.md",
            "modules/Filaments.md",
            "modules/FilamentIO.md",
            "modules/BiotSavart.md",
            "modules/Timestepping.md",
            "modules/Diagnostics.md",
        ],
        "References" => "references.md",
    ],
    warnonly = [:missing_docs],  # TODO can we remove this?
    plugins = [bib],
    repo,
)

##

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
