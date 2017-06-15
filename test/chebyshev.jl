using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

#Chebyshev
facts("Chebyshev") do
for T in (Float32, Float64, Complex64, Complex128)
    context("Matrix{$T}") do
    A=convert(Matrix{T}, randn(n,n))
    T<:Complex && (A+=convert(Matrix{T}, im*randn(n,n)))
    A=A+A'
    A=A'*A #Construct SPD matrix
    b=convert(Vector{T}, randn(n))
    T<:Complex && (b+=convert(Vector{T}, im*randn(n)))
    b=b/norm(b)
    tol = 0.1 #For some reason Chebyshev is very slow
    v = eigvals(A)
    mxv = maximum(v)
    mnv = minimum(v)
    x_cheby, c_cheby= chebyshev(A, b, mxv+(mxv-mnv)/100, mnv-(mxv-mnv)/100, tol=tol, maxiter=10^5, log=true)
    @fact c_cheby.isconverged --> true
    @fact norm(A*x_cheby-b) --> less_than(tol)
    end
end
end
