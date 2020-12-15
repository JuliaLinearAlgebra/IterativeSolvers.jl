module TestMINRES

using IterativeSolvers
using Test
using Random
using SparseArrays
using LinearAlgebra
using LinearMaps

@testset "MINRES" begin

function hermitian_problem(T, n)
    B = rand(T, n, n) + n * I
    A = B + B'
    x = ones(T, n)
    b = B * x
    A, x, b
end

function skew_hermitian_problem(T, n)
    B = rand(T, n, n) + n * I
    A = B - B'
    x = ones(T, n)
    b = A * x
    A, x, b
end

Random.seed!(123)
n = 15


@testset "Hermitian Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A, x, b = hermitian_problem(T, n)
    reltol = sqrt(eps(real(T)))
    x0 = rand(T, n)

    x1, hist1 = minres(A, b, maxiter = 10n, reltol = reltol, log = true)
    x2, hist2 = minres!(x0, A, b, maxiter = 10n, reltol = reltol, log = true)

    @test isa(hist1, ConvergenceHistory)
    @test norm(b - A * x1) / norm(b) ≤ reltol
    @test hist1.isconverged
    @test norm(b - A * x2) / norm(b) ≤ reltol
    @test x2 == x0
end

@testset "Skew-Hermitian Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A, x, b = skew_hermitian_problem(T, n)
    reltol = sqrt(eps(real(T)))
    x_approx, hist = minres(A, b, skew_hermitian = true, maxiter = 10n, reltol = reltol, log = true)

    @test norm(b - A * x_approx) / norm(b) ≤ reltol
    @test hist.isconverged
end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float32, Float64, ComplexF32, ComplexF64), Ti in (Int64, Int32)
    A = let
        B = sprand(n, n, 2 / n)
        B + B' + I
    end

    x = ones(T, n)
    b = A * x
    reltol = sqrt(eps(real(T)))

    x_approx, hist = minres(A, b, maxiter = 10n, reltol = reltol, log = true)

    @test norm(b - A * x_approx) / norm(b) ≤ reltol
    @test hist.isconverged
end

@testset "Termination criterion" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        A = T[ 2 -1  0
              -1  2 -1
               0 -1  2]
        n = size(A, 2)
        b = ones(T, n)
        x0 = A \ b
        perturbation = 10 * sqrt(eps(real(T))) * T[(-1)^i for i in 1:n]

        # If the initial residual is small and a small relative tolerance is used,
        # many iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = minres!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = minres!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end

end

end # module
