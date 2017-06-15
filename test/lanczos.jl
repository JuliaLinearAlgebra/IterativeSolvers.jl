using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

#Lanczos methods

type MyOp{T}
    buf::Matrix{T}
end

import Base: *, size, eltype

*(A::MyOp, x::AbstractVector) = A.buf*x
size(A::MyOp, i) = size(A.buf, i)
eltype(A::MyOp) = eltype(A.buf)

facts("eiglancz") do
for T in (Float32, Float64)
    context("Matrix{$T}") do
    A = convert(Matrix{T}, randn(n,n))
    A = A + A' #Symmetric
    v = eigvals(A)

    eval_lanczos, c_lanczos = eiglancz(A, log=true)
    @fact c_lanczos.isconverged --> true
    @fact norm(v - eval_lanczos) --> less_than(√eps(T))
    end

    context("Op{$T}") do
        A = MyOp(convert(Matrix{T}, randn(5,5)) |> t -> t + t')
        v = eigvals(Symmetric(A.buf))
        eval_lanczos, c_lanczos = eiglancz(A, log=true)
        @fact c_lanczos.isconverged --> true
        @fact norm(v - eval_lanczos) --> less_than(√eps(T))
    end
end
end
