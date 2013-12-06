using IterativeSolvers
using Base.Test
include("getDivGrad.jl")

# small full system
A = [4 1; 1 4]
rhs = [2;2]
x,flag,relres,iter,resvec = cg(A,rhs,1e-15)
@test norm(A*x-rhs)/norm(rhs) <= 1e-15

# CG: test sparse Laplacian
A = getDivGrad(32,32,32)
Af(x) = A*x
L = tril(A)
D = diag(A)
U = triu(A)
JAC(x) = D.\x
SGS(x) = L\(D.*(U\x))

rhs = randn(size(A,1))
tol = 1e-5
# tests with A being matrix
xCG,flagCG,relresCG,iterCG,resvecCG       = cg(A,rhs,tol,100)
xJAC,flagJAC,relresJAC,iterJAC,resvecJAC  = cg(A,rhs,tol,100,JAC)
xSGS,flagSGS,relresSGS,iterSGS,resvecSGS  = cg(A,rhs,tol,100,SGS)
# tests with A being function
xCGmf,flagCG,relresCG,iterCG,resvecCG       = cg(Af,rhs,tol,100)
xJACmf,flagJAC,relresJAC,iterJAC,resvecJAC  = cg(Af,rhs,tol,100,JAC)
xSGSmf,flagSGS,relresSGS,iterSGS,resvecSGS  = cg(Af,rhs,tol,100,SGS)
# tests with random starting guess
xCGr,flagCGr,relresCGr,iterCGr,resvecCGr       = cg(Af,rhs,tol,100,1,randn(size(rhs)))
xJACr,flagJACr,relresJACr,iterJACr,resvecJACr  = cg(Af,rhs,tol,100,JAC,randn(size(rhs)))
xSGSr,flagSGSr,relresSGSr,iterSGSr,resvecSGSr  = cg(Af,rhs,tol,100,SGS,randn(size(rhs)))

# test relative residuals
@test norm(A*xCG-rhs)/norm(rhs) <= tol
@test norm(A*xSGS-rhs)/norm(rhs) <= tol
@test norm(A*xJAC-rhs)/norm(rhs) <= tol
@test norm(A*xCGmf-rhs)/norm(rhs) <= tol
@test norm(A*xSGSmf-rhs)/norm(rhs) <= tol
@test norm(A*xJACmf-rhs)/norm(rhs) <= tol
@test norm(A*xCGr-rhs)/norm(rhs) <= tol
@test norm(A*xSGSr-rhs)/norm(rhs) <= tol
@test norm(A*xJACr-rhs)/norm(rhs) <= tol
# preconditioners should at least not increase number of iter
@test iterJAC==iterCG
@test iterSGS<=iterJAC


