include("getDivGrad.jl")

# small full system
N=10
A = randn(N,N)
A = A'*A
rhs = randn(N)
x, = cg(A,rhs,tol=1e-15)
@test_approx_eq norm(A*x-rhs) N*sqrt(1e-15)

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
xCG, = cg(A,rhs,tol=tol,maxIter=100)
xJAC, = cg(A,rhs,tol=tol,maxIter=100,Preconditioner=JAC)
xSGS, = cg(A,rhs,tol=tol,maxIter=100,Preconditioner=SGS)
# tests with A being function
xCGmf, = cg(Af,rhs,tol=tol,maxIter=100)
xJACmf, = cg(Af,rhs,tol=tol,maxIter=100,Preconditioner=JAC)
xSGSmf, = cg(Af,rhs,tol=tol,maxIter=100,Preconditioner=SGS)
# tests with random starting guess
xCGr, = cg(Af,rhs,tol=tol,maxIter=100,Preconditioner=1,x=randn(size(rhs)))
xJACr, = cg(Af,rhs,tol=tol,maxIter=100,Preconditioner=JAC,x=randn(size(rhs)))
xSGSr, = cg(Af,rhs,tol=tol,maxIter=100,Preconditioner=SGS,x=randn(size(rhs)))

# test relative residuals
@test norm(A*xCG-rhs) <= tol
@test norm(A*xSGS-rhs) <= tol
@test norm(A*xJAC-rhs) <= tol
@test norm(A*xCGmf-rhs) <= tol
@test norm(A*xSGSmf-rhs) <= tol
@test norm(A*xJACmf-rhs) <= tol
@test norm(A*xCGr-rhs) <= tol
@test norm(A*xSGSr-rhs) <= tol
@test norm(A*xJACr-rhs) <= tol
# preconditioners should at least not increase number of iter
@test iterJAC==iterCG
@test iterSGS<=iterJAC

