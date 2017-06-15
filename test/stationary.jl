using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

n = 10
m = 6
srand(1234321)

##################
# Linear solvers #
##################

facts("Stationary solvers") do
    for T in (Float32, Float64, Complex64, Complex128)
        context("$T") do

        A=convert(Matrix{T}, randn(n,n))
        T<:Complex && (A+=convert(Matrix{T}, im*randn(n,n)))
        A+=A'
        A+=n*maximum(abs(A[:]))*diagm(ones(T,n)) #Diagonally dominant
        b=convert(Vector{T}, randn(n))
        T<:Complex && (b+=convert(Vector{T}, im*randn(n)))
        b=b/norm(b)
        x0=convert(Vector{T}, randn(n))
        T<:Complex && (x0+=convert(Vector{T}, im*randn(n)))
        x = A\b
        for solver in [jacobi, gauss_seidel]
            xi, ci=solver(A, b, maxiter=n^4, log=true)
            @fact ci.isconverged --> true
            @fact norm(x-xi) --> less_than(n^3*eps(typeof(real(b[1]))))
        end
        for solver in [jacobi!, gauss_seidel!]
            xi, ci=solver(copy(x0), A, b, maxiter=n^4, log=true)
            @fact ci.isconverged --> true
            @fact norm(x-xi) --> less_than(n^3*eps(typeof(real(b[1]))))
        end

        ω = 0.5
        for solver in [sor, ssor]
            xi, ci=solver(A, b, ω, maxiter=n^4, log=true)
            @fact ci.isconverged --> true
            @fact norm(x-xi) --> less_than(n^3*eps(typeof(real(b[1]))))
        end
        for solver in [sor!, ssor!]
            xi, ci=solver(copy(x0), A, b, ω, maxiter=n^4, log=true)
            @fact ci.isconverged --> true
            @fact norm(x-xi) --> less_than(n^3*eps(typeof(real(b[1]))))
        end

        end
    end
end
