using IterativeSolvers
using Base.Test

function randSPD(T, n)
    A = rand(T, n, n) + n * I
    A' * A
end

function approx_eigenvalue_bounds(A)
    λs = eigvals(A)
    mnv, mxv = minimum(real(λs)), maximum(real(λs))
    Δ = (mxv - mnv) / 100
    mnv - Δ, mxv + Δ
end

#Chebyshev
@testset "Chebyshev" begin

n = 10
srand(1234321)

@testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A = randSPD(T, n)
    b = rand(T, n)
    tol = √(eps(real(T)))
    
    # Without a preconditioner
    begin
        λ_min, λ_max = approx_eigenvalue_bounds(A)
        x0 = rand(n)
        x, history = chebyshev(A, b, λ_min, λ_max, tol=tol, maxiter=10n, log=true)
        @test isa(history, ConvergenceHistory)
        @test history.isconverged
        @test norm(A * x - b) / norm(b) ≤ tol
    end

    # With an initial guess
    begin
        λ_min, λ_max = approx_eigenvalue_bounds(A)
        x0 = rand(T, n)
        x, history = chebyshev!(x0, A, b, λ_min, λ_max, tol=tol, maxiter=10n, log=true)
        @test isa(history, ConvergenceHistory)
        @test history.isconverged
        @test x == x0
        @test norm(A * x - b) / norm(b) ≤ tol
    end

    # With a preconditioner
    begin
        B = randSPD(T, n)
        B_fact = cholfact!(B)
        λ_min, λ_max = approx_eigenvalue_bounds(B_fact \ A)
        x, history = chebyshev(A, b, λ_min, λ_max, Pl = B_fact, tol=tol, maxiter=10n, log=true)
        @test history.isconverged
        @test norm(A * x - b) / norm(b) ≤ tol
    end
end
end
