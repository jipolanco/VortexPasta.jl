using Documenter
using Documenter: Remotes
using DocumenterCitations
using Literate: Literate

using VortexPasta

# We load a few submodules to make @ref's work in Literate-generated files.
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.FilamentIO
using VortexPasta.PredefinedCurves
using VortexPasta.Diagnostics

# For some reason we need to explicitly load these packages here to avoid failures with
# doctests and @example blocks on Gitlab.
using Rotations
using StaticArrays
using SpecialFunctions
using GLMakie
using CairoMakie

# These seem to be needed for GLMakie on Gitlab
using JpegTurbo
using ImageMagick

DocMeta.setdocmeta!(
    VortexPasta.PaddedArrays,
    :DocTestSetup,
    :(using VortexPasta.PaddedArrays),
)

DocMeta.setdocmeta!(
    VortexPasta.PredefinedCurves,
    :DocTestSetup,
    quote
        using StaticArrays: SVector, SMatrix
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

function make_all(; generate_tutorials = true,)
    repo = Gitlab("https://gitlab.in2p3.fr/jipolanco/VortexPasta.jl")
    warnonly = [:missing_docs]  # TODO can we remove this?

    bib = CitationBibliography(
        joinpath(@__DIR__, "src", "biblio.bib");
        style = :authoryear,
    )

    tutorials = String[
        "01-vortex_ring.jl",
        "02-kelvin_waves.jl",
    ]
    if !generate_tutorials
        empty!(tutorials)
        # Some cross-references will be broken if we don't generate the tutorials
        push!(warnonly, :cross_references)
    end

    for name ∈ tutorials
        Literate.markdown(
            joinpath(@__DIR__, "literate", "tutorials", name),
            joinpath(@__DIR__, "src", "tutorials");
            documenter = true,
        )
    end

    Makie.set_theme!()  # reset theme from a possible previous run
    tutorials_gen = map(tutorials) do name
        joinpath("tutorials", replace(name, r"\.jl$" => ".md"))
    end

    makedocs(;
        sitename = "VortexPasta",
        format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true",
            repolink = Remotes.repourl(repo),
            edit_link = "master",
            assets = [
                "assets/fonts.css",
            ],
            size_threshold_ignore = [
                "modules/Filaments.md",  # this page has too much content so it's relatively large
            ],
            mathengine = KaTeX(),
        ),
        modules = [VortexPasta],
        pages = [
            "index.md",
            "Tutorials" => tutorials_gen,
            "Explanations" => [
                "methods/VFM.md",
                "methods/Ewald.md",
            ],
            "Modules" => [
                "modules/PaddedArrays.md",
                "modules/PredefinedCurves.md",
                "modules/CellLists.md",
                "modules/BasicTypes.md",
                "modules/Quadratures.md",
                "modules/Filaments.md",
                "modules/FilamentIO.md",
                "modules/FindNearbySegments.md",
                "modules/BiotSavart.md",
                "modules/Reconnections.md",
                "modules/Timestepping.md",
                "modules/Diagnostics.md",
            ],
            "References" => "references.md",
        ],
        pagesonly = true,
        warnonly,
        plugins = [bib],
        repo,
    )
end

make_all()

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
