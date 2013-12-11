using IterativeSolvers
using Base.Test
n=10
m=6

for T in (Float32, Float64, Complex64, Complex128)
    A=convert(Matrix{T}, randn(n,n))
    T<:Complex && (A+=convert(Matrix{T}, im*randn(n,n)))
    A=A+A' #Symmetric/Hermitian
    v = eigvals(A)
    
    ## Simple methods
    
    #Power iteration
    eval_big = maximum(v) > abs(minimum(v)) ? maximum(v) : minimum(v)
    eval_pow = eigvals_power(A; tol=sqrt(eps(real(one(T)))), maxiter=2000)[1].val
    @test_approx_eq_eps eval_big eval_pow (iseltype(T,Complex)?2:1)*n^2*cond(A)*eps(real(one(T)))
    
    #Inverse iteration
    eval_rand = v[1+int(rand()*(n-1))] #Pick random eigenvalue
    eval_ii = eigvals_ii(A, eval_rand*(1+(rand()-.5)/n); tol=sqrt(eps(real(one(T)))), maxiter=2000)[1].val
    @test_approx_eq_eps eval_rand eval_ii (iseltype(T,Complex)?2:1)*n^2*cond(A)*eps(real(one(T)))
    
    #Rayleigh quotient iteration
    #XXX broken?
    #l = eigvals_rqi(A, eigvals_rand, 2000, sqrt(eps())).val
    #@test_approx_eq eigvals_rand l
end

for T in (Float32, Float64)
    A=convert(Matrix{T}, randn(n,n))
    A=A+A' #Symmetric
    v = eigvals(A)

    #Lanczos methods
    
    #Lanczos eigenvalues computation
    eval_lanczos, = eigvals_lanczos(A)
    @test_approx_eq v eval_lanczos
end


for T in (Float32, Float64)
    B = convert(Matrix{T}, randn(n, m))
    v = svdvals(B)

    #Golub-Kahan-Lanczos singular values computation
    sv_gkl = svdvals_gkl(B)
    @test_approx_eq v sv_gkl
end

#Conjugate gradients
#XXX failing include("cg.jl")

#GMRES
for T in (Float32, Float64, Complex64, Complex128)
    A = convert(Matrix{T}, randn(n,n))
    L = convert(Matrix{T}, randn(n,n))
    R = convert(Matrix{T}, randn(n,n))
    b = convert(Vector{T}, randn(n))
    if T <: Complex
        A += im*convert(Matrix{T}, randn(n,n))
        L += im*convert(Matrix{T}, randn(n,n))
        R += im*convert(Matrix{T}, randn(n,n))
        b += im*convert(Vector{T}, randn(n))
    end
    
    #GMRES
    x_gmres = gmres(A, b, L, R)
    @test_approx_eq A*x_gmres b
end

for T in (Float64, Complex128)
    #GMRES on sparse inputs
    A = sprandn(n,n,0.5)
    L = sprandn(n,n,0.5)
    R = sprandn(n,n,0.5)
    b = convert(Vector{T}, randn(n))
    if T <: Complex
        A += im*sprandn(n,n,0.5)
        L += im*sprandn(n,n,0.5)
        R += im*sprandn(n,n,0.5)
        b += im*randn(n)
    end
    x_gmres = gmres(A, b, L, R)
    @test_approx_eq A*x_gmres b
end

