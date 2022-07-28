module TestQMR

using IterativeSolvers
using Test
using Random
using LinearAlgebra
using SparseArrays

@testset "QMR" begin
# Nb: convergence tests for QMR can be brittle because the residual is not exactly tracked
# during iteration, rather a (tight) estimate

n = 10
m = 6
rng = Random.MersenneTwister(123)

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(rng, T, n, n) + n * I
    b = rand(rng, T, n)
    reltol = √eps(real(T))*10

    x, history = qmr(A, b, log=true)
    @test isa(history, ConvergenceHistory)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(rng, T, n, n, 0.5) + n * I
    b = rand(rng, T, n)
    reltol = √eps(real(T))

    x, history = qmr(A, b, log=true, reltol=reltol)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ 3reltol # TODO: Should maybe not require the 3?
end

@testset "Maximum number of iterations" begin
    x, history = qmr(rand(5, 5), rand(5), log=true, maxiter=2)
    @test history.iters == 2
    @test length(history[:resnorm]) == 2
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
        x, ch = qmr!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = qmr!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end

end

end # module
