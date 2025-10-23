using Documenter
using DocumenterVitepress
using Documenter: Remotes
using DocumenterCitations
using Literate: Literate
using Downloads: Downloads

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
using DSP: DSP

# These seem to be needed for GLMakie on Gitlab
using JpegTurbo
using ImageMagick

DocMeta.setdocmeta!(
    VortexPasta.PaddedArrays,
    :DocTestSetup,
    :(using VortexPasta.PaddedArrays),
)

DocMeta.setdocmeta!(
    VortexPasta.Containers,
    :DocTestSetup,
    :(using VortexPasta.Containers),
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
    VortexPasta.Filaments,
    :DocTestSetup,
    quote
        using StaticArrays: SVector
        using VortexPasta: VortexPasta
        using VortexPasta.Filaments
    end,
)

# Adapted from https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/docs/make.jl
const with_liveserver = "liveserver" ∈ ARGS

if with_liveserver
    using Revise
    Revise.revise()
end

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

# Try to detect repo from environment variables automatically defined by the CI platform.
const REPO = if haskey(ENV, "GITLAB_CI") && haskey(ENV, "CI_PROJECT_URL")  # https://docs.gitlab.com/ee/ci/variables/predefined_variables.html
    Gitlab(ENV["CI_PROJECT_URL"])
else
    # Default: point to github repo
    Remotes.GitHub("jipolanco", "VortexPasta.jl")
end

function make_all(; draft = false,)
    repo = Remotes.repourl(REPO)

    bib = CitationBibliography(
        joinpath(@__DIR__, "src", "biblio.bib");
        style = :authoryear,
    )

    tutorials = String[
        "01-vortex_ring.jl",
        "02-kelvin_waves.jl",
    ]
    if draft
        empty!(tutorials)
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

    # assets = Documenter.HTMLWriter.HTMLAsset[]

    # Try to download latest version of simpleanalytics script.
    try
        script = "public/sa.js"  # see also config.mts where we include this file in <head>
        dst = joinpath(@__DIR__, "src", script)
        Downloads.download("https://scripts.simpleanalyticscdn.com/latest.js", dst)
        # attributes = Dict(:async => "", Symbol("data-collect-dnt") => "true")
        # push!(assets, asset(script; attributes, islocal = true))
    catch e
        @warn "Failed downloading asset" e
    end

    makedocs(;
        sitename = "VortexPasta",
        authors = "Juan Ignacio Polanco",
        format = DocumenterVitepress.MarkdownVitepress(;
            repo,
            devbranch = "master",
            devurl = "dev",
            # collapselevel =  1,
            # assets,
            # size_threshold_ignore = [
            #     "modules/Filaments.md",  # this page has too much content so it's relatively large
            # ],
            # mathengine = KaTeX(),
        ),
        modules = [VortexPasta],
        pages = [
            "index.md",
            "Tutorials" => tutorials_gen,
            "Explanations" => [
                "methods/VFM.md",
                "methods/Ewald.md",
            ],
            "Tips and tricks" => [
                "tips/parallelisation.md",
                "tips/gpu.md",
            ],
            "Modules" => String[
                "modules/PaddedArrays.md",
                "modules/PredefinedCurves.md",
                "modules/CellLists.md",
                "modules/Quadratures.md",
                "modules/Filaments.md",
                "modules/FilamentIO.md",
                "modules/FindNearbySegments.md",
                "modules/Constants.md",
                "modules/BiotSavart.md",
                "modules/SyntheticFields.md",
                "modules/Forcing.md",
                "modules/Containers.md",
                "modules/Reconnections.md",
                "modules/Timestepping.md",
                "modules/Diagnostics.md",
            ],
            "References" => "references.md",
        ],
        draft = draft || with_liveserver,
        pagesonly = true,
        warnonly = true,
        plugins = [bib],
        repo,
        doctest = false,
    )
end

make_all(; draft = false,)

@show readdir(joinpath(@__DIR__, "build"))
@show readdir(joinpath(@__DIR__, "build", "1"))

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
if haskey(ENV, "GITHUB_REPOSITORY")  # if we're on github
    DocumenterVitepress.deploydocs(;
        repo = "github.com/jipolanco/VortexPasta.jl",
        branch = "gh-pages",
        devbranch = "master",
        forcepush = true,
        push_preview = true,
    )
end
