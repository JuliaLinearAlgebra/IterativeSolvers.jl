using IterativeSolvers, Base.Test

tol = 1e-10
N = 10
M = 100
A = sprandn(M, N, 0.2)
b = randn(M)

# sparse matrix
(x, converged) = rkaczmarz(A'A, A'b, tol = tol)
@test norm(A'A*x - A'b) <= tol

# full matrix
(x, converged) = rkaczmarz(full(A'A), A'b, tol = tol)
@test norm(A'A*x - A' * b) <= tol