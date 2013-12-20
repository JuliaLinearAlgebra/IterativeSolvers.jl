require("../src/IterativeSolvers.jl")
using IterativeSolvers
using Base.Test
const n=10
const m=6
srand(1234321)

##################
# Linear solvers #
##################

#Stationary solvers
for T in (Float32, Float64, Complex64, Complex128)
    A=convert(Matrix{T}, randn(n,n))
    T<:Complex && (A+=convert(Matrix{T}, im*randn(n,n)))
    A+=A'
    A+=n*maximum(abs(A[:]))*diagm(ones(T,n)) #Diagonally dominant
    b=convert(Vector{T}, randn(n))
    T<:Complex && (b+=convert(Vector{T}, im*randn(n)))
    x0=convert(Vector{T}, randn(n))
    T<:Complex && (x0+=convert(Vector{T}, im*randn(n)))
    x = A\b
    for guess in {nothing, x0}
        for solver in [jacobi, gauss_seidel]
            xi, ci=solver(A, b, guess, maxiter=n^4)
            @test ci.isconverged
            @test_approx_eq_eps x xi n^3*eps(typeof(real(b[1])))
        end

        ω = 0.5
        for solver in [sor, ssor]
            xi, ci=solver(A, b, ω, guess, maxiter=n^4)
            @test ci.isconverged
            @test_approx_eq_eps x xi n^3*eps(typeof(real(b[1])))
        end
    end
end

#Conjugate gradients
include("cg.jl")

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
    
    x_gmres, c_gmres= gmres(A, b, L, R)
    @test c_gmres.isconverged
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
    x_gmres, c_gmres= gmres(A, b, L, R)
    @test c_gmres.isconverged
    @test_approx_eq A*x_gmres b
end

#Chebyshev
for T in (Float32, Float64, Complex64, Complex128)
    A=convert(Matrix{T}, randn(n,n))
    T<:Complex && (A+=convert(Matrix{T}, im*randn(n,n)))
    A=A+A'
    A=A'*A #Construct SPD matrix
    b=convert(Vector{T}, randn(n))
    T<:Complex && (b+=convert(Vector{T}, im*randn(n)))
    tol = 0.1 #For some reason Chebyshev is very slow
    v = eigvals(A)
    x_cheby, c_cheby= chebyshev(A, b, minimum(v), maximum(v), tol=tol, maxiter=10^5)
    @test c_cheby.isconverged
    @test_approx_eq_eps A*x_cheby b tol
end

#######################
# Eigensystem solvers #
#######################

#Simple eigensolvers
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


#Lanczos methods
    
#Lanczos eigenvalues computation
for T in (Float32, Float64)
    A=convert(Matrix{T}, randn(n,n))
    A=A+A' #Symmetric
    v = eigvals(A)

    eval_lanczos, c_lanczos = eigvals_lanczos(A)
    T==Float64 && @test c_lanczos.isconverged #XXX Lanczos needs to be made more robust for Float32
    @test_approx_eq v eval_lanczos
end

#Golub-Kahan-Lanczos singular values computation
for T in (Float32, Float64)
    B = convert(Matrix{T}, randn(n, m))
    v = svdvals(B)

    sv_gkl = svdvals_gkl(B)
    @test_approx_eq v sv_gkl
end


