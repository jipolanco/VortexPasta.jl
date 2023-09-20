using Documenter
using Documenter: Remotes
using VortexPasta
using Rotations  # loading this here seems to fix doctest issues on Gitlab

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

makedocs(;
    sitename = "VortexPasta",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = Remotes.repourl(repo),
        edit_link = "master",
        size_threshold_ignore = [
            "modules/Filaments.md",  # this page has too much content so it's relatively large
        ],
    ),
    modules = [VortexPasta],
    pages = [
        "index.md",
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
        ]
    ],
    warnonly = [:missing_docs],  # TODO remove
    repo,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
