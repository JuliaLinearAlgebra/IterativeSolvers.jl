using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

facts("MINRES") do

n = 15

for T in (Float32, Float64, Complex64, Complex128)
    context("Matrix{$T}") do
        # Some well-conditioned symmetric matrix for testing
        A = let 
            B = rand(T, n, n) + n * eye(n)
            B + B'
        end

        x = ones(T, n)
        b = A * x
        tol = sqrt(eps(real(T)))

        x_approx, hist = minres(A, b, maxiter = n + 1, tol = tol, log = true)

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

        x_approx, hist = minres(A, b, maxiter = n + 1, tol = tol, log = true)

        @fact norm(b - A * x_approx) / norm(b) --> less_than_or_equal(tol)
        @fact hist.isconverged --> true
    end
end

end