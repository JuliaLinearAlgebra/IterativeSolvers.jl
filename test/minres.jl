using IterativeSolvers
using Base.Test
using LinearMaps

@testset "MINRES" begin

function hermitian_problem(T, n)
    B = rand(T, n, n) + n * I
    A = B + B'
    x = ones(T, n)
    b = B * x
    A, x, b
end

function skew_hermitian_problem(T, n)
    B = rand(T, n, n) + n * I
    A = B - B'
    x = ones(T, n)
    b = A * x
    A, x, b
end

srand(123)
n = 15


@testset "Hermitian Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A, x, b = hermitian_problem(T, n)
    tol = sqrt(eps(real(T)))

    x_approx, hist = minres(A, b, maxiter = 10n, tol = tol, log = true)

    @test isa(hist, ConvergenceHistory)
    @test norm(b - A * x_approx) / norm(b) ≤ tol
    @test hist.isconverged
end

@testset "Skew-Hermitian Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A, x, b = skew_hermitian_problem(T, n)
    tol = sqrt(eps(real(T)))
    x_approx, hist = minres(A, b, skew_hermitian = true, maxiter = 10n, tol = tol, log = true)

    @test norm(b - A * x_approx) / norm(b) ≤ tol
    @test hist.isconverged
end

@testset "SparseMatrixCSC{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A = let
        B = sprand(n, n, 2 / n)
        B + B' + I
    end

    x = ones(T, n)
    b = A * x
    tol = sqrt(eps(real(T)))

    x_approx, hist = minres(A, b, maxiter = 10n, tol = tol, log = true)

    @test norm(b - A * x_approx) / norm(b) ≤ tol
    @test hist.isconverged
end

end