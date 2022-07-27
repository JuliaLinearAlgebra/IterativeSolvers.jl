module TestLALQMR

using IterativeSolvers
using Test
using LinearMaps
using LinearAlgebra
using Random
using SparseArrays

#LALQMR
@testset "LALQMR" begin

Random.seed!(1234321)
n = 10

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(T, n, n)
    b = rand(T, n)
    F = lu(A)
    reltol = √eps(real(T))

    x, history = lalqmr(A, b, log=true, maxiter=10, reltol=reltol);
    @test isa(history, ConvergenceHistory)

    # Left exact preconditioner
    #x, history = lalqmr(A, b, Pl=F, maxiter=1, reltol=reltol, log=true)
    #@test history.isconverged
    #@test norm(F \ (A * x - b)) / norm(b) ≤ reltol

    # Right exact preconditioner
    #x, history = lalqmr(A, b, Pl=Identity(), Pr=F, maxiter=1, reltol=reltol, log=true)
    #@test history.isconverged
    #@test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "SparseMatrixCSC{$T, $Ti}" for T in (Float64, ComplexF64), Ti in (Int64, Int32)
    A = sprand(T, n, n, 0.5) + I
    b = rand(T, n)
    F = lu(A)
    reltol = √eps(real(T))

    x, history = lalqmr(A, b, log = true, maxiter = 10);
    @test norm(A * x - b) / norm(b) ≤ reltol

    # Left exact preconditioner
    #x, history = lalqmr(A, b, Pl=F, maxiter=1, log=true)
    #@test history.isconverged
    #@test norm(F \ (A * x - b)) / norm(b) ≤ reltol

    # Right exact preconditioner
    #x, history = lalqmr(A, b, Pl = Identity(), Pr=F, maxiter=1, log=true)
    #@test history.isconverged
    #@test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "Block Creation {$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    # Guaranteed to create blocks during Lanczos process
    # This satisfies the condition that in the V-W sequence, the first
    # iterates are orthogonal: <Av - v<A, v>, Atv - v<At, v>> under transpose inner product
    dl = fill(one(T), n-1)
    du = fill(one(T), n-1)
    d = fill(one(T), n)
    dl[1] = -1
    A = Tridiagonal(dl, d, du)
    b = fill(zero(T), n)
    b[2] = 1.0

    reltol = √eps(real(T))
    
    x, history = lalqmr(A, b, log = true)
    @test norm(A * x - b) / norm(b) ≤ reltol
end

@testset "Linear operator defined as a function" begin
    A = LinearMap(cumsum!, 100; ismutating=true)
    b = rand(100)
    reltol = 1e-5

    x = lalqmr(A, b; reltol=reltol, maxiter=2000)
    @test norm(A * x - b) / norm(b) ≤ reltol
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
        x, ch = lalqmr!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = lalqmr!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end
end

end # module