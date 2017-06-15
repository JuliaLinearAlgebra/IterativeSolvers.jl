using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

n = 10
m = 6
srand(1234567)

#IDRs
facts("idrs") do

for T in (Float32, Float64, Complex64, Complex128)
    # without smoothing
    context("Matrix{$T}, without residual smoothing") do

    A = convert(Matrix{T}, randn(n,n))
    b = convert(Vector{T}, randn(n))
    if T <: Complex
        A += im*convert(Matrix{T}, randn(n,n))
        b += im*convert(Vector{T}, randn(n))
    end
    b = b/norm(b)

    x_idrs, c_idrs = idrs(A, b, log=true)
    @fact c_idrs.isconverged --> true
    @fact norm(A*x_idrs - b) --> less_than(√eps(real(one(T))))
    end

    # with smoothing
    context("Matrix{$T}, with residual smoothing") do

    A = convert(Matrix{T}, randn(n,n))
    b = convert(Vector{T}, randn(n))
    if T <: Complex
        A += im*convert(Matrix{T}, randn(n,n))
        b += im*convert(Vector{T}, randn(n))
    end
    b = b/norm(b)

    x_idrs, c_idrs = idrs(A, b; smoothing=true, log=true)
    @fact c_idrs.isconverged --> true
    @fact norm(A*x_idrs - b) --> less_than(√eps(real(one(T))))
    end
end

for T in (Float64, Complex128)
    context("SparseMatrixCSC{$T}") do
    A = sprandn(n,n,0.5)+0.001*eye(T,n,n)
    b = convert(Vector{T}, randn(n))
    if T <: Complex
        A += im*sprandn(n,n,0.5)
        b += im*randn(n)
    end
    b = b / norm(b)

    x_idrs, c_idrs= idrs(A, b, log=true)
    @fact c_idrs.isconverged --> true
    @fact norm(A*x_idrs - b) --> less_than(√eps(real(one(T))))
    end
end
end
