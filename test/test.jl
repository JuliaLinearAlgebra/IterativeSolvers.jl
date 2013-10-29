include("../src/IterativeSolvers.jl")
using IterativeSolvers
using Base.Test
n=10
A=randn(n,n)
A=A+A' #Real and symmetric

v = eigvals(A)

println("Power iteration")
k = maximum(v) > abs(minimum(v)) ? maximum(v) : minimum(v)
l = ev_power(A, 2000, sqrt(eps()))[1]
println([k l])
println("Deviation: ", abs(k-l))
@test_approx_eq_eps k l sqrt(eps())

println("Lanczos eigenvalues computation")
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

