module TestGMRES

using IterativeSolvers
using Test
using LinearMaps
using LinearAlgebra
using Random
using SparseArrays

#GMRES
@testset "GMRES" begin

rng = Random.Xoshiro(1234)
n = 10

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(rng, T, n, n) + I
    b = rand(rng, T, n)
    F = lu(A)
    reltol = √eps(real(T))

    # Test optimality condition: residual should be non-increasing
    x, history = gmres(A, b, log=true, restart=3, maxiter=10, reltol=reltol);
    @test isa(history, ConvergenceHistory)
    @test all(diff(history[:resnorm]) .<= 0.0)

    # Left exact preconditioner
    x, history = gmres(A, b, Pl=F, maxiter=1, restart=1, reltol=reltol, log=true)
    @test history.isconverged
    @test norm(F \ (A * x - b)) / norm(b) ≤ reltol

    # Right exact preconditioner
    x, history = gmres(A, b, Pl=Identity(), Pr=F, maxiter=1, restart=1, reltol=reltol, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(rng, T, n, n, 0.5) + I
    b = rand(rng, T, n)
    F = lu(A)
    reltol = √eps(real(T))

    # Test optimality condition: residual should be non-increasing
    x, history = gmres(A, b, log = true, restart = 3, maxiter = 10);
    @test all(diff(history[:resnorm]) .<= 0.0)

    # Left exact preconditioner
    x, history = gmres(A, b, Pl=F, maxiter=1, restart=1, log=true)
    @test history.isconverged
    @test norm(F \ (A * x - b)) / norm(b) ≤ reltol

    # Right exact preconditioner
    x, history = gmres(A, b, Pl = Identity(), Pr=F, maxiter=1, restart=1, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "Linear operator defined as a function" begin
    A = LinearMap(cumsum!, 100; ismutating=true)
    b = rand(rng, 100)
    reltol = 1e-5

    x = gmres(A, b; reltol=reltol, maxiter=2000)
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "Off-diagonal in hessenberg matrix exactly zero" begin
    A = Matrix(1.0I, 2, 2)
    b = [1.0, 2.2]
    x = gmres(A, b)
    @test all(x .== b)
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
        x, ch = gmres!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = gmres!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end
end

end # module
