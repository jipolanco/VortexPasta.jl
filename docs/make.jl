using Documenter
using VortexFilamentEwald

DocMeta.setdocmeta!(
    VortexFilamentEwald.Filaments,
    :DocTestSetup,
    :(using VortexFilamentEwald.Filaments),
)

# doctest(VortexFilamentEwald; fix = true)

makedocs(
    sitename = "VortexFilamentEwald",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [VortexFilamentEwald],
    pages = [
        "index.md",
        "Modules" => [
            "modules/BasicTypes.md",
            "modules/Quadratures.md",
            "modules/Filaments.md",
            "modules/BiotSavart.md",
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
