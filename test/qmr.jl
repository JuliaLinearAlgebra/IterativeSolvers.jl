using IterativeSolvers
using Test
using LinearAlgebra
using SparseArrays
using Random

@testset ("QMR") begin
    Random.seed!(1234321)
    n = 20

    @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
        A = rand(T, n, n) + 15I
        x = ones(T, n)
        b = A * x

        tol = √eps(real(T))

        # Solve without preconditioner
        x1, his1 = qmr(A, b, log = true, tol = tol)
        @test isa(his1, ConvergenceHistory)
        @test norm(A * x1 - b) / norm(b) ≤ tol

        # With an initial guess
        x_guess = rand(T, n)
        x2, his2 = qmr!(x_guess, A, b, log = true, tol = tol)
        @test isa(his2, ConvergenceHistory)
        @test norm(A * x2 - b) / norm(b) ≤ tol

        # Do an exact LU decomp of a nearby matrix
        F = lu(A + rand(T, n, n))
        x3, his3 = qmr(A, b, Pl = F, log = true, tol = tol)
        @test norm(A * x3 - b) / norm(b) ≤ tol
    end
end
