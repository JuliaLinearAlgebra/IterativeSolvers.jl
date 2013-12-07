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
l = ev_power(A, 2000, sqrt(eps()))[1].val
println([k l])
println("Deviation: ", abs(k-l))
@test_approx_eq_eps k l sqrt(eps())

println("Inverse iteration")
ev_rand = v[1+int(rand()*(n-1))] #Pick random eigenvalue
l = ev_ii(A, ev_rand*(1+(rand()-.5)/n), 2000, sqrt(eps()))[1].val
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
w, = eigvals_lanczos(A)
println([v w])
println("Deviation: ", norm(v-w))
#XXX failing! @test_approx_eq v w



m = 6
B = randn(n, m)

println("Golub-Kahan-Lanczos singular values computation")
v = svdvals(B)
w = svdvals_gkl(B)
println([v w])
println("Deviation: ", norm(v-w))
@test_approx_eq v w

#Conjugate gradients
#XXX failing include("cg.jl")

#GMRES
n = 10;
for T in (Float64,Complex{Float64})
    A = rand(T,n,n)
    L = rand(T,n,n)
    R = rand(T,n,n)
    b = rand(T,n)

    println("GMRES $T")
    x = gmres(A, b; M1 = L, M2 = R)
    println([A*x b])
    println("Deviation: ", norm(A*x-b))
    @test_approx_eq A*x b

    println("GMRES Sparse $T")
    A = sparse(A);
    L = sparse(L);
    R = sparse(R);
    x = gmres(A, b; M1 = L, M2 = R)
    println([A*x b])
    println("Deviation: ", norm(A*x-b))
    @test_approx_eq A*x b
end

