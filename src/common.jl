import  Base: eltype, eps, length, ndims, real, size, *,
        A_mul_B!, Ac_mul_B, Ac_mul_B!

export  A_mul_B

#### Type-handling
"""
    Adivtype(A, b)
Determine type of the division of an element of `b` against an element of `A`:
`typeof(one(eltype(b))/one(eltype(A)))`
"""
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))

"""
    Amultype(A, x)
Determine type of the multiplication of an element of `b` with an element of `A`:
`typeof(one(eltype(A))*one(eltype(x)))`
"""
Amultype(A, x) = typeof(one(eltype(A))*one(eltype(x)))
if VERSION < v"0.4.0-dev+6068"
    real{T<:Real}(::Type{Complex{T}}) = T
    real{T<:Real}(::Type{T}) = T
end
eps{T<:Real}(::Type{Complex{T}}) = eps(T)

"""
    randx(A, b)
Build a random unitary vector `Vector{T}`, where `T` is `Adivtype(A,b)`.
"""
function randx(A, b)
    T = Adivtype(A, b)
    x = initrand!(Array(T, size(A, 2)))
end

"""
    zerox(A, b)
Build a zeros vector `Vector{T}`, where `T` is `Adivtype(A,b)`.
"""
function zerox(A, b)
    T = Adivtype(A, b)
    x = zeros(T, size(A, 2))
end

#### Numerics
"""
    initrand!(v)
Overwrite `v` with a random unitary vector of the same length.
"""
function initrand!(v::Vector)
    _randn!(v)
    nv = norm(v)
    for i = 1:length(v)
        v[i] /= nv
    end
    v
end
_randn!(v::Array{Float64}) = randn!(v)
_randn!(v) = copy!(v, randn(length(v)))

#### Errors
export PosSemidefException

type PosSemidefException <: Exception
    msg :: AbstractString
    PosSemidefException(msg::AbstractString="Matrix was not positive semidefinite") = new(msg)
end

#### Support for linear-operators-defined-as-functions
#
# Suppose you have a function implementing multiplication of b by A.
# This function can have any syntax, but for the purposes of
# illustration let's suppose it's defined as
#   mulbyA!(output, b, Adata)
# Where Adata might be some parameters that your function needs.
# You can represent it as a linear operator using
#   A = MatrixFcn{T}(m, n, (output, b) -> mulbyA!(output, b, Adata))
# where T is the "element type" of Adata.
# Note that there are a couple of requirements:
#   - mulbyA! stores the result in the pre-allocated output
#   - mulbyA! should also return output as its sole return value

# If the algorithm also needs multiplication by A', use MatrixCFcn
# instead and initialize it with both functions.

export MatrixFcn, MatrixCFcn

abstract AbstractMatrixFcn{T}

type MatrixFcn{T} <: AbstractMatrixFcn{T}
    m::Int
    n::Int
    mul::Function
end

type MatrixCFcn{T} <: AbstractMatrixFcn{T}
    m::Int
    n::Int
    mul::Function
    mulc::Function
end

MatrixFcn(A::AbstractMatrix) = MatrixFcn{eltype(A)}(size(A, 1), size(A, 2), (output, x)-> A_mul_B!(output, A, x))
MatrixCFcn(A::AbstractMatrix) = MatrixFcn{eltype(A)}(size(A, 1), size(A, 2), (output, x)->A_mul_B!(output, A, x), (output, b) -> Ac_mul_B(output, A, b))

eltype{T}(op::AbstractMatrixFcn{T}) = T

ndims(op::AbstractMatrixFcn) = 2

size(op::AbstractMatrixFcn, dim::Integer) = (dim == 1) ? op.m :
                                            (dim == 2) ? op.n : 1
size(op::AbstractMatrixFcn) = (op.m, op.n)

length(op::AbstractMatrixFcn) = op.m*op.n

*(op::AbstractMatrixFcn, b) = A_mul_B(op, b)

A_mul_B{R,S}(op::AbstractMatrixFcn{R}, b::AbstractVector{S}) = A_mul_B!(Array(promote_type(R,S), op.m), op, b)

A_mul_B!(output, op::AbstractMatrixFcn, b) = op.mul(output, b)

Ac_mul_B{R,S}(op::MatrixCFcn{R}, b::AbstractVector{S}) = Ac_mul_B!(Array(promote_type(R,S), op.n), op, b)

Ac_mul_B!(output, op::AbstractMatrixFcn, b) = op.mulc(output, b)
