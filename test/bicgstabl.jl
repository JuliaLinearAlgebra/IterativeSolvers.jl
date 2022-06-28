module TestBiCGStabl

using IterativeSolvers
using Test
using Random
using LinearAlgebra

@testset ("BiCGStab(l)") begin

Random.seed!(1234321)
n = 20

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    Random.seed!(123) # Issue #316 (test sensitive to the rng)
    A = rand(T, n, n) + 15I
    x = ones(T, n)
    b = A * x

    reltol = √eps(real(T))

    @testset "BiCGStab($l)" for l = (2, 4)
        # Solve without preconditioner
        x1, his1 = bicgstabl(A, b, l, max_mv_products = 100, log = true, reltol = reltol)
        @test isa(his1, ConvergenceHistory)
        @test norm(A * x1 - b) / norm(b) ≤ reltol

        # With an initial guess
        x_guess = rand(T, n)
        x2, his2 = bicgstabl!(x_guess, A, b, l, max_mv_products = 100, log = true, reltol = reltol)
        @test isa(his2, ConvergenceHistory)
        @test x2 == x_guess
        @test norm(A * x2 - b) / norm(b) ≤ reltol

        # The following tests fails CI on Windows and Ubuntu due to a
        # `SingularException(4)`
        if T == Float32 && (Sys.iswindows() || Sys.islinux())
            continue
        end
        # Do an exact LU decomp of a nearby matrix
        F = lu(A + rand(T, n, n))
        x3, his3 = bicgstabl(A, b, Pl = F, l, max_mv_products = 100, log = true, reltol = reltol)
        @test norm(A * x3 - b) / norm(b) ≤ reltol
    end
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
        x, ch = bicgstabl!(x, A, b, log=true)
        @test 1 ≤ niters(ch) ≤ n÷2 # BiCGStab(2) makes more than one RHS evaluations per iteration

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = bicgstabl!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end
end

end # module
