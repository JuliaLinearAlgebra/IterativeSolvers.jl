using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(123)

facts("MINRES") do

function hermitian_problem(T, n)
    B = rand(T, n, n) + n * eye(n)
    A = B + B'
    x = ones(T, n)
    b = B * x
    A, x, b
end

function skew_hermitian_problem(T, n)
    B = rand(T, n, n) + n * eye(n)
    A = B - B'
    x = ones(T, n)
    b = A * x
    A, x, b
end

for T in (Float32, Float64, Complex64, Complex128)
    n = 15

    context("Hermitian Matrix{$T}") do
        A, x, b = hermitian_problem(T, n)
        tol = sqrt(eps(real(T)))

        x_approx, hist = minres(A, b, maxiter = 10n, tol = tol, log = true)

        @fact norm(b - A * x_approx) / norm(b) --> less_than_or_equal(tol)
        @fact hist.isconverged --> true
    end

    context("Skew-Hermitian Matrix{$T}") do
        A, x, b = skew_hermitian_problem(T, n)
        tol = sqrt(eps(real(T)))
        x_approx, hist = minres(A, b, skew_hermitian = true, maxiter = 10n, tol = tol, log = true)

        @fact norm(b - A * x_approx) / norm(b) --> less_than_or_equal(tol)
        @fact hist.isconverged --> true
    end

    context("SparseMatrixCSC{$T}") do
        A = let
            B = sprand(n, n, 2 / n)
            B + B' + speye(n)
        end

        x = ones(T, n)
        b = A * x
        tol = sqrt(eps(real(T)))

        x_approx, hist = minres(A, b, maxiter = 10n, tol = tol, log = true)

        @fact norm(b - A * x_approx) / norm(b) --> less_than_or_equal(tol)
        @fact hist.isconverged --> true
    end
end

end