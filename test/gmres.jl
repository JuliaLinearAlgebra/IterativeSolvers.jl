using IterativeSolvers
using Base.Test
using LinearMaps


#GMRES
@testset "GMRES" begin

srand(1234321)
n = 10

@testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A = rand(T, n, n)
    b = rand(T, n)
    F = lufact(A)
    tol = √eps(real(T))

    # Test optimality condition: residual should be non-increasing
    x, history = gmres(A, b, log = true, restart = 3, maxiter = 10, tol = tol);
    @test all(diff(history[:resnorm]) .<= 0.0)

    # Left exact preconditioner
    x, history = gmres(A, b, Pl=F, maxiter=1, restart=1, tol=tol, log=true)
    @test history.isconverged
    @test norm(F \ (A * x - b)) / norm(b) ≤ tol

    # Right exact preconditioner
    x, history = gmres(A, b, Pl=Identity(), Pr=F, maxiter=1, restart=1, tol=tol, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol
end

@testset "SparseMatrixCSC{$T}" for T in (Float64, Complex128)
    A = sprand(T, n, n, 0.5) + speye(T, n, n)
    b = rand(T, n)
    F = lufact(A)
    tol = √eps(real(T))

    # Test optimality condition: residual should be non-increasing
    x, history = gmres(A, b, log = true, restart = 3, maxiter = 10);
    @test all(diff(history[:resnorm]) .<= 0.0)

    # Left exact preconditioner
    x, history = gmres(A, b, Pl=F, maxiter=1, restart=1, log=true)
    @test history.isconverged
    @test norm(F \ (A * x - b)) / norm(b) ≤ tol

    # Right exact preconditioner
    x, history = gmres(A, b, Pl = Identity(), Pr=F, maxiter=1, restart=1, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol
end

@testset "Linear operator defined as a function" begin
    A = LinearMap(cumsum!, 100; ismutating=true)
    b = rand(100)
    tol = 1e-5

    x = gmres(A, b; tol=tol, maxiter=2000)
    @test norm(A * x - b) / norm(b) ≤ tol
end
end
