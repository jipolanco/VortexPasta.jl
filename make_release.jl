# To make a new release:
#
# 1. Update project version in Project.toml
# 2. Run this script from a terminal (`julia make_release.jl`)
#
# This script expects LocalRegistry to be installed in the global environment.

using Pkg: Pkg
using LocalRegistry: LocalRegistry

Pkg.activate(@__DIR__)
project = Pkg.project()
@assert project.version !== nothing

# If using TagBot (on github), we let TagBot create the tag, instead of doing it ourselves.
using_tagbot = true

version = string(project.version)  # e.g. "0.7.1"
tag = "v" * version  # e.g. "v0.7.1"

function run_stdout(cmd)
    io = IOBuffer()
    run(cmd, devnull, io)
    String(take!(io))
end

if using_tagbot
    run(`git push`)  # only push latest commits; don't create a new tag
else
    # Check if tag already exists
    tag_exists = !isempty(run_stdout(`git tag -l $tag`))
    if tag_exists
        # Check that the tag points to the current commit.
        hash_tag = run_stdout(`git rev-parse $tag`)
        hash_current = run_stdout(`git rev-parse HEAD`)
        if hash_tag == hash_current
            @info "Tag $tag already exists. Not re-tagging."
        else
            @error "Tag $tag already exists but it points to a different commit!"
        end
    else
        run(`git tag $tag`)
    end
    run(`git push -o ci.skip`)  # push latest commits but don't trigger CI
    run(`git push --tags`)      # push tags (also triggers CI)
end

LocalRegistry.register()
