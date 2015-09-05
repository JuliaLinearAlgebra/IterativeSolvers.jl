using IterativeSolvers, Base.Test

tol = 1e-10
N=10
M = 100
A = sprandn(M,N, 0.2)
x = randn(N)
rhs = A*x

# sparse matrix
(x, converged) = rkaczmarz(A',rhs, tol = tol)
@test norm(A*x - rhs) <= tol

# full matrix
(x, converged) = rkaczmarz(full(A'),rhs, tol = tol)
@test norm(A*x - rhs) <= tol