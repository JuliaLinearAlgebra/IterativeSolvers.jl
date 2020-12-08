using IterativeSolvers
using Test
using Random
using LinearAlgebra
using SparseArrays

@testset "QMR" begin

n = 10
m = 6
Random.seed!(1234567)

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(T, n, n) + n * I
    b = rand(T, n)
    reltol = √eps(real(T))*10

    x, history = qmr(A, b, log=true)
    @test isa(history, ConvergenceHistory)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(T, n, n, 0.5) + n * I
    b = rand(T, n)
    reltol = √eps(real(T))

    x, history = qmr(A, b, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol
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
        perturbation = T[(-1)^i for i in 1:n]

        # If the initial residual is small and a small relative tolerance is used,
        # many iterations are necessary
        x = x0 + sqrt(eps(real(T))) * perturbation
        initial_residual = norm(A * x - b)
        x, ch = qmr!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + 10*sqrt(eps(real(T))) * perturbation
        initial_residual = norm(A * x - b)
        x, ch = qmr!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end

end
