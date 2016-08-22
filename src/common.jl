import  Base: eltype, eps, length, ndims, real, size, *, \, ctranspose,
        A_mul_B!, Ac_mul_B, Ac_mul_B!

export  A_mul_B

# Improve readability of iterative methods
\(f::Function, b) = f(b)
*(f::Function, b) = f(b)

#### Type-handling
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))
Amultype(A, x) = typeof(one(eltype(A))*one(eltype(x)))
if VERSION < v"0.4.0-dev+6068"
    real{T<:Real}(::Type{Complex{T}}) = T
    real{T<:Real}(::Type{T}) = T
end
eps{T<:Real}(::Type{Complex{T}}) = eps(T)

function randx(A, b)
    T = Adivtype(A, b)
    x = initrand!(Array(T, size(A, 2)))
end

function zerox(A, b)
    T = Adivtype(A, b)
    x = zeros(T, size(A, 2))
end

#### Numerics
function update!(x, α::Number, p::AbstractVector)
    for i = 1:length(x)
        x[i] += α*p[i]
    end
    x
end

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

export FuncMatrix

"""
    FuncMatrix{T,A,B}
Represent functions as a matrix.
**Fields**
* `m::Int` = number of columns.
* `n::Int` = number of rows.
* `mul::Function` = `A*b` implementation.
* `cmul::Function` = `A'*b` implementation.
**Constructors**
    FuncMatrix(A)
    FuncMatrix(
        m::Int, n::Int; typ::Type=Float64, ctrans::Bool=false,
        mul::Function=identity, cmul=identity
    )
**Arguments**
* `A::AbstractMatrix` = matrix.
* `m::Int` = number of columns.
* `n::Int` = number of rows.
* `typ::Type = Float64` = `eltype(::FuncMatrix)`.
* `mul::Function = identity` = `A*b` implementation.
* `cmul::Function = identity` = `A'*b` implementation.
"""
type FuncMatrix{T,A,B}
    m::Int
    n::Int
    mul::Function
    cmul::Function
end
function FuncMatrix(
        m::Int, n::Int; typ::Type=Float64,
        mul::Function=identity, cmul::Function=identity
        )
    FuncMatrix{typ,mul!=identity,cmul!=identity}(m, n, mul, cmul)
end
FuncMatrix(A::AbstractMatrix) = FuncMatrix(
    size(A, 1),
    size(A, 2),
    typ = eltype(A),
    mul = (output, x)->A_mul_B!(output, A, x),
    cmul = (output, b) -> Ac_mul_B(output, A, b)
    )

eltype{T}(::FuncMatrix{T}) = T

ndims(::FuncMatrix) = 2

size(op::FuncMatrix) = (op.m, op.n)
size(op::FuncMatrix, dim::Integer) = (dim == 1) ? op.m : (dim == 2) ? op.n : 1

length(op::FuncMatrix) = op.m*op.n

ctranspose{T,A,B}(op::FuncMatrix{T,A,B}) = FuncMatrix{T,B,A}(op.n, op.m, op.cmul, op.mul)

*(op::FuncMatrix, b) = A_mul_B(op, b)

function A_mul_B{R,S}(op::FuncMatrix{R}, b::AbstractVector{S})
    A_mul_B!(Array(promote_type(R,S), op.m), op, b)
end
function A_mul_B{R,S}(op::FuncMatrix{R}, b::AbstractMatrix{S})
    A_mul_B!(Array(promote_type(R,S), op.m, size(b,2)), op, b)
end

A_mul_B!{T}(output, op::FuncMatrix{T,false}, b) = error("A*b not defined")
function A_mul_B!{T}(output, op::FuncMatrix{T,true}, b::AbstractVector)
    op.mul(output, b)
end
function A_mul_B!{T}(output, op::FuncMatrix{T,true}, b::AbstractMatrix)
    columns = [op.mul(output, b[:,i]) for i in 1:size(b,2)]
    hcat(columns...)
end

function Ac_mul_B{R,S}(op::FuncMatrix{R}, b::AbstractVector{S})
    Ac_mul_B!(Array(promote_type(R,S), op.n), op, b)
end
function Ac_mul_B{R,S}(op::FuncMatrix{R}, b::AbstractMatrix{S})
    Ac_mul_B!(Array(promote_type(R,S), op.n, size(b,2)), op, b)
end

Ac_mul_B!{T,A}(output, op::FuncMatrix{T,A,false}, b) = error("A'*b not defined")
function Ac_mul_B!{T,A}(output, op::FuncMatrix{T,A,true}, b::AbstractVector)
    op.cmul(output, b)
end
function Ac_mul_B!{T,A}(output, op::FuncMatrix{T,A,true}, b::AbstractMatrix)
    columns = [op.cmul(output, b[:,i]) for i in 1:size(b,2)]
    hcat(columns...)
end
