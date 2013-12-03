using IterativeSolvers
using Base.Test
include("getDivGrad.jl")

# small full system
A = [4 1; 1 4]; rhs = [2;2];
x,flag,relres,iter,resvec = pcg(A,rhs,1e-15)
@test norm(A*x-rhs)/norm(rhs) <= 1e-15

# CG: test sparse Laplacian
A = getDivGrad(32,32,32)
L = tril(A);
D = diag(A);
U = triu(A);
JAC(x) = D.\x;
SGS(x) = L\(D.*(U\x));

rhs = randn(size(A,1))
tol = 1e-5
tic(); xCG,flagCG,relresCG,iterCG,resvecCG       = pcg(A,rhs,tol,100); toc();
tic(); xJAC,flagJAC,relresJAC,iterJAC,resvecJAC  = pcg(A,rhs,tol,100,JAC); toc();
tic(); xSGS,flagSGS,relresSGS,iterSGS,resvecSGS  = pcg(A,rhs,tol,100,SGS); toc();
@test norm(A*xCG-rhs)/norm(rhs) <= tol
@test norm(A*xSGS-rhs)/norm(rhs) <= tol
@test norm(A*xJAC-rhs)/norm(rhs) <= tol
@test iterJAC==iterCG
@test iterSGS<=iterJAC
