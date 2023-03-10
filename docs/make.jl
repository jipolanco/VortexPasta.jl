using Documenter
using VortexFilamentEwald

DocMeta.setdocmeta!(
    VortexFilamentEwald.Filaments,
    :DocTestSetup,
    :(using VortexFilamentEwald.Filaments),
)

makedocs(
    sitename = "VortexFilamentEwald",
    format = Documenter.HTML(),
    modules = [VortexFilamentEwald],
    pages = [
        "index.md",
        "Modules" => [
            "modules/Filaments.md",
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
