using IterativeSolvers
using Test
using LinearMaps
using LinearAlgebra
using Random
using SparseArrays

#GMRES
@testset "GMRES" begin

Random.seed!(1234321)
n = 10

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(T, n, n)
    b = rand(T, n)
    F = lu(A)
    tol = √eps(real(T))

    # Test optimality condition: residual should be non-increasing
    x, history = gmres(A, b, log = true, restart = 3, maxiter = 10, tol = tol);
    @test isa(history, ConvergenceHistory)
    @test all(diff(history[:resnorm]) .<= 0.0)

    # Left exact preconditioner
    x, history = gmres(A, b, Pl=F, maxiter=1, restart=1, tol=tol, log=true)
    @test history.isconverged
    @test norm(F \ (A * x - b)) / norm(b) ≤ tol

    # Right exact preconditioner
    x, history = gmres(A, b, Pl=Identity(), Pr=F, maxiter=1, restart=1, tol=tol, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol
end

@testset "SparseMatrixCSC{$T}" for T in (Float64, ComplexF64)
    A = sprand(T, n, n, 0.5) + I
    b = rand(T, n)
    F = lu(A)
    tol = √eps(real(T))

    # Test optimality condition: residual should be non-increasing
    x, history = gmres(A, b, log = true, restart = 3, maxiter = 10);
    @test all(diff(history[:resnorm]) .<= 0.0)

    # Left exact preconditioner
    x, history = gmres(A, b, Pl=F, maxiter=1, restart=1, log=true)
    @test history.isconverged
    @test norm(F \ (A * x - b)) / norm(b) ≤ tol

    # Right exact preconditioner
    x, history = gmres(A, b, Pl = Identity(), Pr=F, maxiter=1, restart=1, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol
end

@testset "Linear operator defined as a function" begin
    A = LinearMap(cumsum!, 100; ismutating=true)
    b = rand(100)
    tol = 1e-5

    x = gmres(A, b; tol=tol, maxiter=2000)
    @test norm(A * x - b) / norm(b) ≤ tol
end
end
