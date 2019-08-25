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
    tol = √eps(real(T))*10

    x, history = qmr(A, b, log=true)
    @test isa(history, ConvergenceHistory)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol


end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(T, n, n, 0.5) + n * I
    b = rand(T, n)
    tol = √eps(real(T))

    x, history = qmr(A, b, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol
end

@testset "Maximum number of iterations" begin
    x, history = qmr(rand(5, 5), rand(5), log=true, maxiter=2)
    @test history.iters == 2
    @test length(history[:resnorm]) == 2
end

end
