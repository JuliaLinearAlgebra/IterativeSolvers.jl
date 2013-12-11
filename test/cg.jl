include("getDivGrad.jl")

# small full system
N=10
A = randn(N,N)
A = A'*A
rhs = randn(N)
x, = cg(A,rhs;tol=1e-15)
@test_approx_eq_eps A*x rhs cond(A)*sqrt(1e-15)

# CG: test sparse Laplacian
A = getDivGrad(32,32,32)
Af(x) = A*x
L = tril(A)
D = diag(A)
U = triu(A)
JAC(x) = D.\x
SGS(x) = L\(D.*(U\x))

rhs = randn(size(A,2))
tol = 1e-5
# tests with A being matrix
xCG, = cg(A,rhs;tol=tol,maxiter=100)
xJAC, = cg(A,rhs,JAC;tol=tol,maxiter=100)
xSGS, = cg(A,rhs,SGS;tol=tol,maxiter=100)
# tests with A being function
xCGmf, = cg(Af,rhs;tol=tol,maxiter=100)
xJACmf, = cg(Af,rhs,JAC;tol=tol,maxiter=100)
xSGSmf, = cg(Af,rhs,SGS;tol=tol,maxiter=100)
# tests with random starting guess
xCGr, hCGr = cg(Af,rhs,x->x,randn(size(rhs));tol=tol,maxiter=100)
xJACr, hJACr = cg(Af,rhs,JAC,randn(size(rhs));tol=tol,maxiter=100)
xSGSr, hSGSr = cg(Af,rhs,SGS,randn(size(rhs));tol=tol,maxiter=100)

# test relative residuals
@test_approx_eq_eps A*xCG rhs tol
@test_approx_eq_eps A*xSGS rhs tol
@test_approx_eq_eps A*xJAC rhs tol
@test_approx_eq_eps A*xCGmf rhs tol
@test_approx_eq_eps A*xSGSmf rhs tol
@test_approx_eq_eps A*xJACmf rhs tol
@test_approx_eq_eps A*xCGr rhs tol
@test_approx_eq_eps A*xSGSr rhs tol
@test_approx_eq_eps A*xJACr rhs tol
# preconditioners should at least not increase number of iter
iterCG = length(hCGr.residuals)
iterJAC = length(hJACr.residuals)
iterSGS = length(hSGSr.residuals)
@test iterJAC==iterCG
@test iterSGS<=iterJAC

