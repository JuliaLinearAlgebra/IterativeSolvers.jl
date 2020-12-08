using IterativeSolvers
using Test
using Random
using LinearAlgebra

@testset ("BiCGStab(l)") begin

Random.seed!(1234321)
n = 20

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(T, n, n) + 15I
    x = ones(T, n)
    b = A * x

    tol = √eps(real(T))

    @testset "BiCGStab($l)" for l = (2, 4)
        # Solve without preconditioner
        x1, his1 = bicgstabl(A, b, l, max_mv_products = 100, log = true, tol = tol)
        @test isa(his1, ConvergenceHistory)
        @test norm(A * x1 - b) / norm(b) ≤ tol

        # With an initial guess
        x_guess = rand(T, n)
        x2, his2 = bicgstabl!(x_guess, A, b, l, max_mv_products = 100, log = true, tol = tol)
        @test isa(his2, ConvergenceHistory)
        @test x2 == x_guess
        @test norm(A * x2 - b) / norm(b) ≤ tol

        # The following tests fails CI on Windows and Ubuntu due to a
        # `SingularException(4)`
        if T == Float32 && (Sys.iswindows() || Sys.islinux())
            continue
        end
        # Do an exact LU decomp of a nearby matrix
        F = lu(A + rand(T, n, n))
        x3, his3 = bicgstabl(A, b, Pl = F, l, max_mv_products = 100, log = true, tol = tol)
        @test norm(A * x3 - b) / norm(b) ≤ tol
    end
end
end
