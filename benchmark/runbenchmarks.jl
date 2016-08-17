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
function runBG(bg::BenchmarkGroup, name::AbstractString)
    println("Running...")
    result = minimum(run(bg; verbose=true))
    y = convert(Array{AbstractString,1},collect(keys(result)))
    x = convert(Array{Float64,1},map(y -> time(y),collect(values(result))))
    unit = getunit(maximum(x))
    x = map(y->prettytime(y,unit), x)
    println()
    println(barplot(y,x,title = "$name ($unit)"))
end

# Executes all the suites inside: "Dense", "Sparse" and "Function"
#
#   suite       ->  Suite
function runSuite(suite)
    for name in ["Dense", "Sparse", "Function"]
        println("$name Suite:\n")
        for key in keys(suite[name])
            runBG(suite[name][key],"$name $key")
        end
    end
end

###Run Benchmarks###
println("\nStarting matrix test suites ...\n")
runSuite(Benchmarks.suite)
