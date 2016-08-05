include(joinpath(Pkg.dir("IterativeSolvers"), "benchmark", "BenchmarkUtils.jl"))
using .BenchmarkUtils
include(joinpath(Pkg.dir("IterativeSolvers"), "benchmark", "Benchmarks.jl"))
using .Benchmarks

using BenchmarkTools
using UnicodePlots

# Runs and plots the given benchmarkgroup.
#
#   bg          ->  Suite
#   name        ->  Plot title
function getPlot(bg::BenchmarkGroup, name::AbstractString)
    result = minimum(run(bg))
    y = convert(Array{AbstractString,1},collect(keys(result)))
    x = convert(Array{Float64,1},map(y -> time(y),collect(values(result))))
    unit = getunit(maximum(x))
    x = map(y->prettytime(y,unit), x)
    barplot(y,x,title = "$name ($unit)")
end

# Executes all the suites inside: "Dense", "Sparse" and "Function"
#
#   suite       ->  Suite
function getPlots(suite)
    matrices = Dict()
    for name in ["Dense", "Sparse", "Function"]
        matrices[name] = Dict()
        for key in keys(suite[name])
            matrices[name][key] = getPlot(suite[name][key],"$key")
        end
    end
    matrices
end

# Removes the second bar in unicode plots, because the width of the space change
# in markdown.
function remove_end(uniplot::AbstractString)
    count = 0
    res = ""
    for i in uniplot
        (i == '│') && (count = (count + 1) % 3)
        (count == 2) && (count=0; continue)
        res = res * "$i"
    end
    res
end

matrices = getPlots(Benchmarks.suite)

file = open(joinpath(Pkg.dir("IterativeSolvers"), "docs", "src", "benchmarks.md"),"w")
println(file, "# Benchmarks\n\n")

println(file,
    "This document provides benchmarks of the iterative methods over " *
    "well-known matrices. The matrices are all of size 10. For more benchmark
    options go to the [benchmark](https://github.com/JuliaLang/IterativeSolvers.jl/tree/master/benchmark)
    folder in the package.\n"
    )

println(file, """```@contents
Pages = ["benchmarks.md"]
```""")
println(file)

for stage in ["Dense", "Sparse", "Function"]
    println(file, "## $stage\n")
    for i in keys(matrices[stage])
        plot = matrices[stage][i]
        println(file,
            "### $i\n" *
            "```\n" *
            remove_end("$plot") *
            "```\n"
            )
    end
end

close(file)
