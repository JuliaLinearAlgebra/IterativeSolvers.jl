using IterativeSolvers
using LinearMaps
using Base.Test

include("advection_diffusion.jl")

@testset ("BiCGStab(l)") begin

srand(1234321)
n = 20

@testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A = rand(T, n, n) + 15 * eye(T, n)
    x = ones(T, n)
    b = A * x

    tol = √eps(real(T))

    @testset "BiCGStab($l)" for l = (2, 4)
        # Solve without preconditioner
        x1, his1 = bicgstabl(A, b, l, max_mv_products = 100, log = true, tol = tol)
        @test norm(A * x1 - b) / norm(b) ≤ tol

        # Do an exact LU decomp of a nearby matrix
        F = lufact(A + rand(T, n, n))
        x2, his2 = bicgstabl(A, b, Pl = F, l, max_mv_products = 100, log = true, tol = tol)
        @test norm(A * x2 - b) / norm(b) ≤ tol
    end
end
end
