#!/usr/bin/env julia
#
# Benchmarking script for singular value problems

try
    Pkg.installed("Benchmarks")
catch
    Pkg.clone("https://github.com/johnmyleswhite/Benchmarks.jl")
end

try
    Pkg.installed("PROPACK")
catch
    Pkg.clone("https://github.com/andreasnoack/PROPACK.jl")
end

using Benchmarks
using IterativeSolvers
using JLD
using MAT
using PROPACK

BASEDIR = "florida"

if !isdir(BASEDIR)
    info("Setting up $BASEDIR collection")
    include("setup-florida.jl")
end

"""
Run benchmark function `runbenchmark` against all .mat files in all
subdirectories of BASEDIR
"""
function benchmarkall(BASEDIR)
    for group in readdir(BASEDIR)
        isdir(joinpath(BASEDIR, group)) || continue
        for matrix in readdir(joinpath(BASEDIR, group))
            filename = joinpath(BASEDIR, group, matrix)

            endswith(matrix, ".mat") || continue

        #If we already saved benchmarks, don't rerun
        benchmarkfilename = joinpath(BASEDIR, group, matrix[1:end-3]*"jld")
        isfile(benchmarkfilename) && continue

        #To debug this script, it's useful to run it on small matrices only
        #filesize(filename) < 100000 || continue

        runbenchmark(filename, benchmarkfilename)

        #debug: just try one
        #@goto stop
        end
    end
    @label stop
end


"""
Run several different computations for singular values
"""
function runbenchmark(filename, benchmarkfilename)
    mf = matread(filename)
    if !haskey(mf, "Problem")
        warn("Skipping unknown file $filename: No 'Problem' struct")
        return
    end

    pr = mf["Problem"]
    if !haskey(pr, "A")
        warn("Skipping unknown file $filename: 'Problem' struct has no matrix 'A'")
        return
    end
    A = pr["A"]
    #eltype(A) <: Real || warn("Skipping matrix with unsupported element type $(eltype(A))")
    info(filename*", matrix of type $(typeof(A)) and size $(size(A))")
    m, n = size(A)
    #Choose the same normalized unit vector to start with
    q = randn(n)
    eltype(A) <: Complex && (q += im*randn(n))
    scale!(q, inv(norm(q)))
    qm = randn(m)
    eltype(A) <: Complex && (q += im*randn(m))
    scale!(qm, inv(norm(qm)))

    #Number of singular values to request
    nv = min(m, n, 10)

    #Maximum number of iterations
    maxiter = max(m, n)

    #Tolerance, however the algorithm chooses to interpret this
    #Set convergence criterion to sqrt eps
    tol = √eps(real(one(eltype(A))))

    info("svds (eigs on [0 A; A' 0])")
    #See #12890 :(
    #b_svds = @benchmark svds(A, v0=q, nsv=nv, tol=tol, maxiter=maxiter)
    b_svds = @benchmark svds(A, nsv=nv, tol=tol, maxiter=maxiter)

    info("eigs on A'A or AA'")
    makeata(A) = (size(A, 1)≥ size(A,2) ? A*A' : A'A)
    B = makeata(A)
    b_ata = @benchmark makeata(A)
    b_eigs = if m≥n
        B = A'A #Size n x n
        b_eigs = @benchmark eigs(B, v0=q, nev=nv, tol=tol, maxiter=maxiter)
    else
        B = A*A' #Size m x m
        #eigs cannot use the same starting vector as everyone else, needs R^m not R^n
        b_eigs = @benchmark eigs(B, v0=qm, nev=nv, tol=tol, maxiter=maxiter)
    end

    info("tsvd (PROPACK)")
    b_tsvd = try
        if m≥n #PROPACK implicitly assumes the input is tall
            @benchmark tsvdvals(A, initvec=qm, kmax=maxiter, k=nv, tolin=tol)
        else
            @benchmark tsvdvals(A', initvec=q, kmax=maxiter, k=nv, tolin=tol)
        end
    catch exc
        println("Exception: $exc")
    end

    info("GKL with thick restart using Ritz values")
    b_tr = try
        @benchmark svdl(A, nv, v0=q, tol=tol, reltol=tol)
    catch exc
        println("Exception: $exc")
    end

    info("GKL with thick restart using harmonic Ritz values")
        b_trh = try
            @benchmark svdl(A, nv, v0=q, tol=tol, reltol=tol, method=:harmonic)
        catch exc
            println("Exception: $exc")
        end

    data = Dict(
        "svds" => b_svds,
        "ata"  => b_ata,
        "eigs" => b_eigs,
        "tsvd" => b_tsvd,
        "tr"   => b_tr,
        "trh"  => b_trh
    )
    println("Timings:")
    for (k, v) in data
        println(k)
        display(v)
        println()
    end
    JLD.save(benchmarkfilename, data, compress=true)
end

benchmarkall(BASEDIR)
