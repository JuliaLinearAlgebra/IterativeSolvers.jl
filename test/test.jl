include("../src/IterativeSolvers.jl")
using IterativeSolvers
using Base.Test
n=10
A=randn(n,n)
A=A+A' #Real and symmetric

println("Lanczos eigenvalues computation")
v = eigvals(A)
w = eigvals_lanczos(A)
println([v w])
println("Deviation: ", norm(v-w))
@test_approx_eq v w

m = 6
B = randn(n, m)

println("Golub-Kahan-Lanczos singular values computation")
v = svdvals(B)
w = svdvals_gkl(B)
println([v w])
println("Deviation: ", norm(v-w))
@test_approx_eq v w

