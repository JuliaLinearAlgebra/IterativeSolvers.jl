include("../src/IterativeSolvers.jl")
using IterativeSolvers
using Base.Test
n=10
A=randn(n,n)
A=A+A' #Real and symmetric

v = eigvals(A)

#Simple methods
println("Power iteration")
k = maximum(v) > abs(minimum(v)) ? maximum(v) : minimum(v)
l = ev_power(A, 2000, sqrt(eps())).val
println([k l])
println("Deviation: ", abs(k-l))
@test_approx_eq_eps k l sqrt(eps())

println("Inverse iteration")
ev_rand = v[1+int(rand()*(n-1))] #Pick random eigenvalue
l = ev_ii(A, ev_rand*(1+(rand()-.5)/n), 2000, sqrt(eps())).val
println([ev_rand l])
println("Deviation: ", abs(ev_rand-l))
@test_approx_eq_eps ev_rand l sqrt(eps())

#XXX broken?
#println("Rayleigh quotient iteration")
#l = ev_rqi(A, ev_rand, 2000, sqrt(eps())).val
#println([ev_rand l])
#println("Deviation: ", abs(ev_rand-l))
#@test_approx_eq_eps ev_rand l sqrt(eps())



#Lanczos methods
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

