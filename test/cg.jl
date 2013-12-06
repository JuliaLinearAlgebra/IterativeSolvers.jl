using IterativeSolvers
using Base.Test
include("getDivGrad.jl")

# small full system
A = [4 1; 1 4]
rhs = [2;2]
x, = pcg(A,rhs,1e-15)
@test norm(A*x-rhs) <= 1e-15

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
xCG, = pcg(A,rhs,tol,100)
xJAC, = pcg(A,rhs,tol,100,JAC)
xSGS, = pcg(A,rhs,tol,100,SGS)
# tests with A being function
xCGmf, = pcg(Af,rhs,tol,100)
xJACmf, = pcg(Af,rhs,tol,100,JAC)
xSGSmf, = pcg(Af,rhs,tol,100,SGS)
# tests with random starting guess
xCGr, = pcg(Af,rhs,tol,100,1,randn(size(rhs)))
xJACr, = pcg(Af,rhs,tol,100,JAC,randn(size(rhs)))
xSGSr, = pcg(Af,rhs,tol,100,SGS,randn(size(rhs)))

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

