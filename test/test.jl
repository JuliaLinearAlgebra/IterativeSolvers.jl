include("../src/IterativeSolvers.jl")
using IterativeSolvers
using Base.Test
n=10
A=randn(n,n)
A=A+A' #Real and symmetric

v = eigvals(A)
w = eigvals_lanczos(A)
println([v w])
println("Deviation: ", norm(eigvals(A) - eigvals_lanczos(A)))
@test_approx_eq v w
