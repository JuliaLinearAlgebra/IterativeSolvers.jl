module TestSimpleEigensolver

using IterativeSolvers
using Test
using LinearMaps
using LinearAlgebra
using Random

@testset "Simple Eigensolvers" begin

rng = Random.Xoshiro(1234)
n = 10

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)

    A = rand(rng, T, n, n) + I
    A = A' * A
    λs = eigvals(A)

    tol = (eltype(T) <: Complex ? 2 : 1) * n^2 * cond(A) * eps(real(one(T)))

    ## Simple methods

    @testset "Power iteration" begin
        λ, x, history = powm(A; tol=tol, maxiter=10n, log=true)
        @test isa(history, ConvergenceHistory)
        @test λs[end] ≈ λ
        @test norm(A * x - λ * x) ≤ tol
    end

    @testset "Inverse iteration" begin
        # Set a target near the middle eigenvalue
        idx = div(n, 2)
        σ = T(0.75 * λs[idx] + 0.25 * λs[idx + 1])

        # Construct F = inv(A - σI) "matrix free"
        # Make sure we use complex arithmetic everywhere,
        # because of the follow bug in base: https://github.com/JuliaLang/julia/issues/22683
        F = lu(complex(A) - UniformScaling(σ))
        Fmap = LinearMap{complex(T)}((y, x) -> ldiv!(y, F, x), size(A, 1), ismutating = true)

        λ, x, history = invpowm(Fmap; shift=σ, tol=tol, maxiter=10n, log=true)

        @test isa(history, ConvergenceHistory)
        @test norm(A * x - λ * x) ≤ tol
        @test λ ≈ λs[idx]

    end
end
end

end # module
