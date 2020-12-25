module TestStationary

import LinearAlgebra.SingularException
import IterativeSolvers: DiagonalIndices, FastLowerTriangular,
       FastUpperTriangular, forward_sub!, backward_sub!, OffDiagonal,
       gauss_seidel_multiply!, StrictlyUpperTriangular, StrictlyLowerTriangular

using IterativeSolvers
using Test
using Random
using LinearAlgebra
using SparseArrays

@testset "Stationary solvers" begin

n = 10
m = 6
ω = 1.2
Random.seed!(1234322)

@testset "SparseMatrix{$T, $Ti}" for T in (Float32, Float64, ComplexF32, ComplexF64), Ti in (Int64, Int32)
    @testset "Sparse? $sparse" for sparse = (true, false)
        # Diagonally dominant
        if sparse
            A = sprand(T, n, n, 4 / n) + 2n * I
        else
            A = rand(T, n, n) + 2n * I
        end

        b = rand(T, n)
        x0 = rand(T, n)
        tol = √eps(real(T))

        for solver in (jacobi, gauss_seidel)
            xi = solver(A, b, maxiter=2n)
            @test norm(b - A * xi) / norm(b) ≤ tol
        end

        for solver in (jacobi!, gauss_seidel!)
            xi = solver(copy(x0), A, b, maxiter=2n)
            @test norm(b - A * xi) / norm(b) ≤ tol
        end

        for solver in (sor, ssor)
            xi = solver(A, b, ω, maxiter=2n)
            @test norm(b - A * xi) / norm(b) ≤ tol
        end

        for solver in (sor!, ssor!)
            xi = solver(copy(x0), A, b, ω, maxiter=2n)
            @test norm(b - A * xi) / norm(b) ≤ tol
        end
    end
end

@testset "SOR and Gauss-Seidel" begin
    A = sprand(10, 10, 4/10) + 4I
    x = ones(10)
    b = A * x

    # Gauss-Seidel and SOR should coincide when ω = 1
    gauss_seidel_it = IterativeSolvers.gauss_seidel_iterable(zeros(10), A, b, maxiter = 5)
    sor_it = IterativeSolvers.sor_iterable(zeros(10), A, b, 1.0, maxiter = 5)

    for _ = zip(gauss_seidel_it, sor_it)
        @test gauss_seidel_it.x ≈ sor_it.x
    end
end

@testset "Should throw for singular diagonals" begin
    A = [0.0 1.0; 1.0 0.0]
    B = sparse(A)
    b = rand(2)

    for matrix in (A, B)
        for solver in (jacobi, gauss_seidel)
            @test_throws LinearAlgebra.SingularException solver(matrix, b)
        end

        for solver in (sor, ssor)
            @test_throws LinearAlgebra.SingularException solver(A, b, 0.5)
        end
    end
end

@testset "Diagonal" begin
    A = sprand(10, 10, .3) + 10I
    B = sparse([0. 0. 0.; 0. 1. 0.; 0. 0. 1.])
    D = DiagonalIndices(A)
    @test A.nzval[D.diag] == Diagonal(A).diag
    @test_throws SingularException DiagonalIndices(B)

    x = ones(10)
    y = zeros(10)
    ldiv!(y, D, x)
    @test y ≈ Diagonal(A) \ ones(10)
end

@testset "Forward substitution" begin
    A = sprand(10, 10, .3) + 10I
    D = DiagonalIndices(A)
    L = FastLowerTriangular(A, D)

    x = rand(10)
    y = copy(x)

    forward_sub!(L, x)

    @test x ≈ LowerTriangular(A) \ y
end

@testset "Forward substitution with update" begin
    B = rand(3, 3) + 10I
    A = sparse(B)
    D = DiagonalIndices(A)
    L = FastLowerTriangular(A, D)
    x = rand(3)
    y = ones(3)
    z = copy(x)
    α = 2.0
    β = 3.0

    forward_sub!(α, L, x, β, y)

    z[1] = β * y[1] + α * z[1] / B[1,1]
    z[2] = β * y[2] + α * (z[2] - B[2,1] * z[1]) / B[2, 2]
    z[3] = β * y[3] + α * (z[3] - B[3,1] * z[1] - B[3,2] * z[2]) / B[3, 3]

    @test x ≈ z
end

@testset "Backward substitution" begin
    A = sprand(10, 10, .3) + 10I
    D = DiagonalIndices(A)
    L = FastUpperTriangular(A, D)

    x = rand(10)
    y = copy(x)

    backward_sub!(L, x)

    @test x ≈ UpperTriangular(A) \ y
end

@testset "Backward substitution with update" begin
    B = rand(3, 3) + 10I
    A = sparse(B)
    D = DiagonalIndices(A)
    L = FastUpperTriangular(A, D)
    x = rand(3)
    y = ones(3)
    z = copy(x)
    α = 2.0
    β = 3.0

    backward_sub!(α, L, x, β, y)

    z[3] = β * y[3] + α * z[3] / B[3,3]
    z[2] = β * y[2] + α * (z[2] - B[2,3] * z[3]) / B[2, 2]
    z[1] = β * y[1] + α * (z[1] - B[1,2] * z[2] - B[1,3] * z[3]) / B[1, 1]

    @test x ≈ z
end

@testset "mul! with OffDiagonal type" begin
    A = sprand(10, 10, .3) + 2I
    D = DiagonalIndices(A)
    O = OffDiagonal(A, D)

    x = rand(10)
    y = ones(10)
    mul!(1.0, O, x, 0.0, y)
    @test y ≈ A * x - Diagonal(A) * x

    x = rand(10)
    y = ones(10)
    mul!(1.0, O, x, 1.0, y)
    @test y ≈ A * x - Diagonal(A) * x + ones(10)

    x = rand(10)
    y = ones(10)
    mul!(2.0, O, x, 3.0, y)
    @test y ≈ 2 * (A * x - Diagonal(A) * x) + 3 * ones(10)
end

@testset "Gauss-Seidel multiplication upper" begin
    A = sprand(10, 10, .3) + 2I
    D = DiagonalIndices(A)
    U = StrictlyUpperTriangular(A, D)

    x = rand(10)
    y = copy(x)
    b = rand(10)

    # x ← b - U * x
    gauss_seidel_multiply!(-1.0, U, x, 1.0, b, x)

    @test x ≈ b - triu(A, 1) * y
end

@testset "Gauss-Seidel multiplication lower" begin
    A = sprand(10, 10, .3) + 2I
    D = DiagonalIndices(A)
    L = StrictlyLowerTriangular(A, D)

    x = rand(10)
    y = copy(x)
    b = rand(10)

    # x ← b - L * x
    gauss_seidel_multiply!(-1.0, L, x, 1.0, b, x)

    @test x ≈ b - tril(A, -1) * y
end
end

end # module
