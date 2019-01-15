using IterativeSolvers
using LinearMaps
using Test
using LinearAlgebra
using SparseArrays
using Random

import LinearAlgebra.ldiv!

include("laplace_matrix.jl")

struct JacobiPrec
    diagonal
end

ldiv!(y, P::JacobiPrec, x) = y .= x ./ P.diagonal

@testset "Conjugate Gradients" begin

Random.seed!(1234321)

@testset "Small full system" begin
    n = 10

    @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
        A = rand(T, n, n)
        A = A' * A + I
        b = rand(T, n)
        tol = √eps(real(T))

        r = cg(A, b; tol=tol, maxiter=2n, log=true)
        x, ch = r.x, r.history
        @test isa(ch, ConvergenceHistory)
        @test norm(A*x - b) / norm(b) ≤ tol
        @test ch.isconverged

        # If you start from the exact solution, you should converge immediately
        r = cg!(A \ b, A, b; tol=10tol, log=true)
        x, ch = r.x, r.history
        @test niters(ch) ≤ 1
        @test nprods(ch) ≤ 2

        # Test with cholfact should converge immediately
        F = cholesky(A, Val(false))
        r = cg(A, b; Pl=F, log=true)
        x, ch = r.x, r.history
        @test niters(ch) ≤ 2
        @test nprods(ch) ≤ 2

        # All-zeros rhs should give all-zeros lhs
        r = cg(A, zeros(T, n))
        x0 = r.x
        @test x0 == zeros(T, n)
    end
end

@testset "Sparse Laplacian" begin
    A = laplace_matrix(Float64, 10, 2)
    P = JacobiPrec(diag(A))

    rhs = randn(size(A, 2))
    rmul!(rhs, inv(norm(rhs)))
    tol = 1e-5

    @testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, Float32), Ti in (Int64, Int32)
        r = cg(A, rhs; tol=tol, maxiter=100)
        xCG = r.x
        r = cg(A, rhs; Pl=P, tol=tol, maxiter=100)
        xJAC = r.x
        @test norm(A * xCG - rhs) ≤ tol
        @test norm(A * xJAC - rhs) ≤ tol
    end

    Af = LinearMap(A)
    @testset "Function" begin
        r = cg(Af, rhs; tol=tol, maxiter=100)
        xCG = r.x
        r = cg(Af, rhs; Pl=P, tol=tol, maxiter=100)
        xJAC = r.x
        @test norm(A * xCG - rhs) ≤ tol
        @test norm(A * xJAC - rhs) ≤ tol
    end

    @testset "Function with specified starting guess" begin
        x0 = randn(size(rhs))
        r = cg!(copy(x0), Af, rhs; tol=tol, maxiter=100, log=true)
        xCG, hCG = r.x, r.history
        r = cg!(copy(x0), Af, rhs; Pl=P, tol=tol, maxiter=100, log=true)
        xJAC, hJAC = r.x, r.history
        @test norm(A * xCG - rhs) ≤ tol
        @test norm(A * xJAC - rhs) ≤ tol
        @test niters(hJAC) == niters(hCG)
    end
end

@testset "CG with a view" begin
    A = rand(10, 10)
    A = A + A' + 100I
    x = view(rand(10, 2), :, 1)
    b = rand(10)
    r = cg!(x, A, b, log = true)
    x, hist = r.x, r.history
    @test hist.isconverged
end

end
