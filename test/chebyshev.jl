using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

function randSPD(T, n)
    A = rand(T, n, n) + n * eye(T, n)
    A = A' + A
    A = A' * A
end

#Chebyshev
facts("Chebyshev") do
n = 10
for T in (Float32, Float64, Complex64, Complex128)
    context("Matrix{$T}") do
        A = randSPD(T, n)
        b = rand(T, n)
        b /= norm(b)

        tol = sqrt(eps(real(T)))
        mnv, mxv = eigmin(A), eigmax(A)
        Δ = (mxv - mnv) / 100

        x, history = chebyshev(A, b, mnv - Δ, mxv + Δ, tol=tol, maxiter=10n, log=true)
        @fact history.isconverged --> true
        @fact norm(A * x - b) --> less_than(tol)

        context("Preconditioned") do
            B = randSPD(T, n)
            B_fact = cholfact!(B)
            BA = B_fact \ A
            λs = eigvals(BA)
            mnv, mxv = minimum(real(λs)), maximum(real(λs))
            Δ = (mxv - mnv) / 100

            x, history = chebyshev(A, b, mnv - Δ, mxv + Δ, Pl = B_fact, tol=tol, maxiter=10n, log=true)
            @fact history.isconverged --> true
            @fact norm(A * x - b) --> less_than(tol)
        end
    end
end
end
