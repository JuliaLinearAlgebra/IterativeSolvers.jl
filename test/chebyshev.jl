using IterativeSolvers
using Test
using Random
using LinearAlgebra

function randSPD(T, n)
    A = rand(T, n, n) + n * I
    A' * A
end

function approx_eigenvalue_bounds(A)
    λs = eigvals(A)
    mnv, mxv = minimum(real(λs)), maximum(real(λs))
    Δ = (mxv - mnv) / 100
    mnv - Δ, mxv + Δ
end

#Chebyshev
@testset "Chebyshev" begin

n = 10
Random.seed!(1234321)

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = randSPD(T, n)
    b = rand(T, n)
    reltol = √(eps(real(T)))
    abstol = reltol

    @testset "Without preconditioner" begin
        λ_min, λ_max = approx_eigenvalue_bounds(A)
        x0 = rand(n)
        x, history = chebyshev(A, b, λ_min, λ_max, reltol=reltol, maxiter=10n, log=true)
        @test isa(history, ConvergenceHistory)
        @test history.isconverged
        @test norm(A * x - b) / norm(b) ≤ reltol
    end

    @testset "With an initial guess" begin
        λ_min, λ_max = approx_eigenvalue_bounds(A)
        x0 = rand(T, n)
        initial_residual = norm(A * x0 - b)
        x, history = chebyshev!(x0, A, b, λ_min, λ_max, reltol=reltol, maxiter=10n, log=true)
        @test isa(history, ConvergenceHistory)
        @test history.isconverged
        @test x == x0
        @test norm(A * x - b) ≤ reltol * initial_residual

        x0 = rand(T, n)
        x, history = chebyshev!(x0, A, b, λ_min, λ_max, abstol=abstol, reltol=zero(real(T)), maxiter=10n, log=true)
        @test isa(history, ConvergenceHistory)
        @test history.isconverged
        @test x == x0
        @test norm(A * x - b) ≤ 2*abstol
    end

    @testset "With a preconditioner" begin
        B = randSPD(T, n)
        B_fact = cholesky!(B, Val(false))
        λ_min, λ_max = approx_eigenvalue_bounds(B_fact \ A)
        x, history = chebyshev(A, b, λ_min, λ_max, Pl = B_fact, reltol=reltol, maxiter=10n, log=true)
        @test history.isconverged
        @test norm(A * x - b) / norm(b) ≤ reltol
    end
end

@testset "Termination criterion" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        A = T[ 2 -1  0
              -1  2 -1
               0 -1  2]
        n = size(A, 2)
        b = ones(T, n)
        x0 = A \ b
        perturbation = T[(-1)^i for i in 1:n]
        λ_min, λ_max = approx_eigenvalue_bounds(A)

        # If the initial residual is small and a small relative tolerance is used,
        # many iterations are necessary
        x = x0 + sqrt(eps(real(T))) * perturbation
        initial_residual = norm(A * x - b)
        x, ch = chebyshev!(x, A, b, λ_min, λ_max, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + 10*sqrt(eps(real(T))) * perturbation
        initial_residual = norm(A * x - b)
        x, ch = chebyshev!(x, A, b, λ_min, λ_max, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end
end
