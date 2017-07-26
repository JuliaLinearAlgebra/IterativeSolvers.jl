using IterativeSolvers
using LinearMaps
using Base.Test

srand(1234321)
include("poisson_matrix.jl")

@testset "Conjugate Gradients" begin

@testset "Small full system" begin
    n = 10

    @testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
        A = rand(T, n, n)
        A = A' * A + I
        b = rand(T, n)
        tol = √eps(real(T))

        x,ch = cg(A, b; tol=tol, maxiter=2n, log=true)
        @test norm(A*x - b) / norm(b) ≤ tol
        @test ch.isconverged

        # If you start from the exact solution, you should converge immediately
        x,ch = cg!(A \ b, A, b; tol=10tol, log=true)
        @test niters(ch) ≤ 1
        @test nprods(ch) ≤ 2

        # Test with cholfact should converge immediately
        F = cholfact(A)
        x,ch = cg(A, b; Pl=F, log=true)
        @test niters(ch) ≤ 2
        @test nprods(ch) ≤ 2

        # All-zeros rhs should give all-zeros lhs
        x0 = cg(A, zeros(T, n))
        @test x0 == zeros(T, n)
    end
end

@testset "Sparse Laplacian" begin
    A = poisson_matrix(Float64, 10, 3)
    L = tril(A)
    D = diag(A)
    U = triu(A)

    JAC(x) = D .\ x
    SGS(x) = L \ (D .* (U \ x))

    rhs = randn(size(A, 2))
    rhs /= norm(rhs)
    tol = 1e-5

    @testset "Matrix" begin
        xCG = cg(A, rhs; tol=tol, maxiter=100)
        xJAC = cg(A, rhs; Pl=JAC, tol=tol, maxiter=100)
        xSGS = cg(A, rhs; Pl=SGS, tol=tol, maxiter=100)
        @test norm(A * xCG - rhs) ≤ tol
        @test norm(A * xSGS - rhs) ≤ tol
        @test norm(A * xJAC - rhs) ≤ tol
    end

    Af = LinearMap(A)
    @testset "Function" begin
        xCG = cg(Af, rhs; tol=tol, maxiter=100)
        xJAC = cg(Af, rhs; Pl=JAC, tol=tol, maxiter=100)
        xSGS = cg(Af, rhs; Pl=SGS, tol=tol, maxiter=100)
        @test norm(A * xCG - rhs) ≤ tol
        @test norm(A * xSGS - rhs) ≤ tol
        @test norm(A * xJAC - rhs) ≤ tol
    end

    @testset "Function with specified starting guess" begin
        tol = 1e-4
        x0 = randn(size(rhs))
        xCG, hCG = cg!(copy(x0), Af, rhs; tol=tol, maxiter=100, log=true)
        xJAC, hJAC = cg!(copy(x0), Af, rhs; Pl=JAC, tol=tol, maxiter=100, log=true)
        xSGS, hSGS = cg!(copy(x0), Af, rhs; Pl=SGS, tol=tol, maxiter=100, log=true)
        @test norm(A * xCG - rhs) ≤ tol
        @test norm(A * xSGS - rhs) ≤ tol
        @test norm(A * xJAC - rhs) ≤ tol
        @test niters(hJAC) == niters(hCG)
        @test niters(hSGS) ≤ niters(hJAC)
    end
end

end
