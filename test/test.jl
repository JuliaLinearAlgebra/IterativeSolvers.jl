include("../src/IterativeSolvers.jl")
using IterativeSolvers

n=10
A=randn(n,n)
A=A+A' #Real and symmetric

evals = [eigvals(A) eigvals_lanczos(A, n, true)]
println(evals)

