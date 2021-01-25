module TestCOCG

using IterativeSolvers
using Test
using Random
using LinearAlgebra

@testset ("Conjugate Orthogonal Conjugate Gradient") begin

Random.seed!(1234321)
n = 20

@testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A = rand(T, n, n)
    A = A + transpose(A) + 15I
    x = ones(T, n)
    b = A * x

    reltol = √eps(real(T))

    # Solve without preconditioner
    x1, his1 = cocg(A, b, reltol = reltol, maxiter = 100, log = true)
    @test isa(his1, ConvergenceHistory)
    @test norm(A * x1 - b) / norm(b) ≤ reltol

    # With an initial guess
    x_guess = rand(T, n)
    x2, his2 = cocg!(x_guess, A, b, reltol = reltol, maxiter = 100, log = true)
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
    x3, his3 = cocg(A, b, Pl = F, maxiter = 100, reltol = reltol, log = true)
    @test norm(A * x3 - b) / norm(b) ≤ reltol
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
        x, ch = cocg!(x, A, b, log=true)
        @test 2 ≤ niters(ch) ≤ n

        # If the initial residual is small and a large absolute tolerance is used,
        # no iterations are necessary
        x = x0 + perturbation
        initial_residual = norm(A * x - b)
        x, ch = cocg!(x, A, b, abstol=2*initial_residual, reltol=zero(real(T)), log=true)
        @test niters(ch) == 0
    end
end
end

end # module
