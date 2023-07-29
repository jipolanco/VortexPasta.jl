using Documenter
using VortexPasta

DocMeta.setdocmeta!(
    VortexPasta.PaddedArrays,
    :DocTestSetup,
    :(using VortexPasta.PaddedArrays),
)

DocMeta.setdocmeta!(
    VortexPasta.PredefinedCurves,
    :DocTestSetup,
    :(using VortexPasta.PredefinedCurves),
)

DocMeta.setdocmeta!(
    VortexPasta.BasicTypes,
    :DocTestSetup,
    :(using VortexPasta.BasicTypes),
)

DocMeta.setdocmeta!(
    VortexPasta.Filaments,
    :DocTestSetup,
    :(using VortexPasta.Filaments),
)

# doctest(VortexPasta; fix = true)

makedocs(
    sitename = "VortexPasta",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
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
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
