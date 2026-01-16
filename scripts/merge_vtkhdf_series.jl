const usage = """
Merge multiple *.vtkhdf.series JSON files, typically resulting from consecutive simulations.

Usage:
    
    julia merge_vtkhdf_series.jl [OPTIONS] [file1.vtkhdf.series file2.vtkhdf.series ...]

Options:

    -o    Output file (stdout by default)

In practice, one generally has a list of vtkhdf.series files in the same directory (e.g.
`simulation_run{1,2,3,…}.vtkhdf.series`) for different runs, in which case one can simply do:

    julia merge_vtkhdf_series.jl -o simulation.vtkhdf.series simulation_run*.vtkhdf.series
"""

if "--help" in ARGS
    println(usage)
    exit(0)
end

# Activate temporary environment where packages will be installed.
using Pkg
Pkg.activate(; temp = true)
Pkg.add(["JSON", "StructUtils"])

# See test/filament_io.jl for some details.
using JSON: JSON
using StructUtils: @tags  # for JSON reading

# This is for parsing JSON files using JSON.jl.
struct ParaViewTimeSeriesItem
    name::String
    time::Float64
end

#  https://juliaio.github.io/JSON.jl/dev/reading/#Field-customization-through-tags
@tags struct ParaViewTimeSeriesFile
    file_series_version::String & (json = (name = "file-series-version",),)
    files::Vector{ParaViewTimeSeriesItem}
end

# Sort files according to simulation time.
function sort_files!(tsf::ParaViewTimeSeriesFile)
    (; files) = tsf
    sort!(files; by = f -> f.time)
    tsf
end

function parse_args(args)
    i = firstindex(args) - 1
    output = nothing  # this is the default
    inputs = String[]
    while i < lastindex(args)
        i += 1
        arg = args[i]
        if arg == "-o"
            i < lastindex(args) || error("found -o option, but without a value")
            i += 1
            output = args[i]
        else
            push!(inputs, arg)  # assume this is an input file
        end
    end
    (; output, inputs)
end

function main(args)
    (; output, inputs) = parse_args(args)
    expected_version = "1.0"
    tsf_out = ParaViewTimeSeriesFile(expected_version, ParaViewTimeSeriesItem[])
    for filename in inputs
        endswith(".series")(filename) || error("expected a .series file; got '$filename'")
        @info "Reading from $filename"
        tsf = JSON.parsefile(filename, ParaViewTimeSeriesFile)::ParaViewTimeSeriesFile
        (; file_series_version, files) = tsf
        file_series_version == expected_version || error("expected file-series-version = \"$expected_version\"; found \"$file_series_version\" ($filename)")
        append!(tsf_out.files, files)
    end
    sort_files!(tsf_out)
    if output === nothing
        @info "Writing to stdout"
        JSON.json(stdout, tsf_out; pretty = true)
    else
        @info "Writing to $output"
        JSON.json(output, tsf_out; pretty = true)
    end
    nothing
end

@main
