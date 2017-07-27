using IterativeSolvers
using Base.Test

@testset "Stationary solvers" begin

import Base.LinAlg.SingularException

n = 10
m = 6
ω = 0.5
srand(1234321)

@testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    # Diagonally dominant matrix.
    A = rand(T, n, n) + 2n * I
    b = rand(T, n)
    x0 = rand(T, n)
    x = A \ b
    tol = √eps(real(T))

    for solver in (jacobi, gauss_seidel)
        xi, history = solver(A, b, maxiter=100, tol=tol, log=true)
        @test history.isconverged
        @test norm(b - A * xi) / norm(b) ≤ tol
    end

    for solver in (jacobi!, gauss_seidel!)
        xi, history = solver(copy(x0), A, b, maxiter=100, tol=tol, log=true)
        @test history.isconverged
        @test norm(b - A * xi) / norm(b) ≤ tol
    end

    for solver in (sor, ssor)
        xi, history = solver(A, b, ω, maxiter=100, tol=tol, log=true)
        @test history.isconverged
        @test norm(b - A * xi) / norm(b) ≤ tol
    end

    for solver in (sor!, ssor!)
        xi, history = solver(copy(x0), A, b, ω, maxiter=100, tol=tol, log=true)
        @test history.isconverged
        @test norm(b - A * xi) / norm(b) ≤ tol
    end

    # Check whether the methods throw when the diagonal has zeros
    A = [0.0 1.0; 1.0 0.0]
    b = rand(3, 3)

    for solver in (jacobi, gauss_seidel)
        @fact_throws Base.LinAlg.SingularException solver(A, b)
    end

    for solver in (sor, ssor)
        @fact_throws Base.LinAlg.SingularException solver(A, b, 0.5)
    end
end
end
