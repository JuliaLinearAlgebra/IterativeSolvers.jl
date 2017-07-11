using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

#######################
# Eigensystem solvers #
#######################

facts("simple eigensolvers") do
for T in (Float32, Float64, Complex64, Complex128)
    context("Matrix{$T}") do
    A=convert(Matrix{T}, randn(n,n))
    T<:Complex && (A+=convert(Matrix{T}, im*randn(n,n)))
    A=A+A' #Symmetric/Hermitian

    tol = (eltype(T) <: Complex ?2:1)*n^2*cond(A)*eps(real(one(T)))
    v = eigvals(A)

    ## Simple methods

    context("Power iteration") do
    eval_big = maximum(v) > abs(minimum(v)) ? maximum(v) : minimum(v)
    eval_pow = powm(A; tol=sqrt(eps(real(one(T)))), maxiter=2000)[1]
    @fact norm(eval_big-eval_pow) --> less_than(tol)
    end

    context("Inverse iteration") do
    irnd = ceil(Int, rand()*(n-2))
    eval_rand = v[1+irnd] #Pick random eigenvalue
    # Perturb the eigenvalue by < 1/4 of the distance to the nearest eigenvalue
    eval_diff = min(abs(v[irnd]-eval_rand), abs(v[irnd+2]-eval_rand))
    σ = eval_rand + eval_diff/2*(rand()-.5)
    eval_ii = invpowm(A; shift=σ, tol=sqrt(eps(real(one(T)))), maxiter=2000)[1]
    @fact norm(eval_rand-eval_ii) --> less_than(tol)
    end

    end
end
end
