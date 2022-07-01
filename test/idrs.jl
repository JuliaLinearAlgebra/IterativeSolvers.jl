module TestIDRs

using IterativeSolvers
using Test
using Random
using LinearAlgebra
using SparseArrays

@testset "IDR(s)" begin

n = 10
m = 6
Random.seed!(1234567)

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(T, n, n) + n * I
    b = rand(T, n)
    reltol = √eps(real(T))

    @testset "Without residual smoothing" begin
        x, history = idrs(A, b, reltol=reltol, log=true)
        @test isa(history, ConvergenceHistory)
        @test history.isconverged
        @test norm(A * x - b) / norm(b) ≤ reltol
    end

    # with smoothing
    @testset "With residual smoothing" begin
        x, history = idrs(A, b; reltol=reltol, smoothing=true, log=true)
        @test history.isconverged
        @test norm(A*x - b) / norm(b) ≤ 2reltol # TODO: Should maybe not require the 2?
    end
end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(T, n, n, 0.5) + n * I
    b = rand(T, n)
    reltol = √eps(real(T))

    x, history = idrs(A, b, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "SparseMatrixCSC{$T, $Ti} with preconditioner" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(T, 1000, 1000, 0.1) + 30 * I
    b = rand(T, 1000)
    reltol = √eps(real(T))

    x, history = idrs(A, b, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol

    Apre = lu(droptol!(copy(A), 0.1)) # inexact preconditioner
    xpre, historypre = idrs(A, b, Pl = Apre, log=true)
    @test historypre.isconverged
    @test norm(A * xpre - b) / norm(b) ≤ reltol

    @test isapprox(x, xpre, rtol = 1e-3)
    @test historypre.iters < 0.5history.iters

end

@testset "Maximum number of iterations" begin
    x, history = idrs(rand(5, 5), rand(5), log=true, maxiter=2)
    @test history.iters == 2
    @test length(history[:resnorm]) == 2
end

@testset "Near solution (#222)" begin
    x = rand(5)
    A = rand(5, 5)
    b = rand(5)

    x, history = idrs!(x, A, b, log=true)
    x_new = copy(x)
    x_new, history = idrs!(x_new, A, b, log=true)

    @test x_new ≈ x
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
        x, ch = idrs!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = idrs!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end

end

end # module
