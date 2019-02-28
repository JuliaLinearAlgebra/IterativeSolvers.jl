using IterativeSolvers
using Test
using LinearAlgebra
using SparseArrays
using Random

@testset "QMR" begin

function matlab_example(T,n)
    A = spdiagm(-1 => fill(-1,n-1), 1 => fill(4,n-1))
    b = sum(A,dims=2)
    M1 = spdiagm(-1 => fill(-1/2,n-1), 0 => ones(n))
    M2 = spdiagm(0 => fill(4,n), 1 => fill(-1,n-1))
    x = ones(T,100)
    A,x,b
end

Random.seed!(123)

@testset "MATLAB Example Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A, x, b = matlab_example(T, n)
    tol = 1e-8

    x0 = rand(T,n)

    x1, hist1 = qmr(A, b, maxiter = 15, tol = tol, log = true)
    x2, hist2 = qmr!(x0, A, b, maxiter = 15, tol = tol, log = true)

    @test isa(hist1, ConvergenceHistory)
    @test norm(b - A * x1) / norm(b) ≤ tol
    @test hist1.isconverged
    @test norm(b - A * x2) / norm(b) ≤ tol
    @test x2 == x0
end

end
