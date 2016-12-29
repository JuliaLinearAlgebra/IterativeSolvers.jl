import  Base: eltype, eps, length, ndims, real, size, *, \, ctranspose,
        A_mul_B!, Ac_mul_B, Ac_mul_B!

export  A_mul_B

# Improve readability of iterative methods
\(f::Function, b) = f(b)
*(f::Function, b) = f(b)

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

export LinearMap

"""
    LinearMap{T,A,B}

Represent a finite dimensional linear map as a function describing
its application (a matrix-vector product) and also, optionally, its transpose.

`T` is the element type of the vector produced by the map.

**Fields**

* `m::Int` = number of columns.
* `n::Int` = number of rows.
* `mul::Function` = `A*b` implementation.
* `cmul::Function` = `A'*b` implementation.

**Constructors**

    LinearMap(A)
    LinearMap(m, n)
    LinearMap(typ, m, n)

**Arguments**

* `A::AbstractMatrix` = matrix.
* `m::Int` = number of columns.
* `n::Int` = number of rows.
* `typ::Type` = element type of the vector produced by the map, if not given
`T` defaults to Float64.
* `mul::Function = identity` = `A*b` implementation.
* `cmul::Function = identity` = `A'*b` implementation.

"""
#A is true if mul is different from identity, otherwise it is false.
#The same happens for B with cmul.
type LinearMap{T,A,B}
    m::Int
    n::Int
    mul::Function
    cmul::Function
end
function LinearMap(
        m::Int, n::Int;
        mul::Function=identity, cmul::Function=identity
        )
    LinearMap{Float64,mul!=identity,cmul!=identity}(m, n, mul, cmul)
end
function LinearMap(
        typ::Type, m::Int, n::Int;
        mul::Function=identity, cmul::Function=identity
        )
    LinearMap{typ,mul!=identity,cmul!=identity}(m, n, mul, cmul)
end
LinearMap(A::AbstractMatrix) = LinearMap(
    eltype(A),
    size(A, 1),
    size(A, 2),
    mul = (output, x)->A_mul_B!(output, A, x),
    cmul = (output, b) -> Ac_mul_B(output, A, b)
    )

eltype{T}(::LinearMap{T}) = T

ndims(::LinearMap) = 2

size(op::LinearMap) = (op.m, op.n)
size(op::LinearMap, dim::Integer) = (dim == 1) ? op.m : (dim == 2) ? op.n : 1

length(op::LinearMap) = op.m*op.n

ctranspose{T,A,B}(op::LinearMap{T,A,B}) = LinearMap{T,B,A}(op.n, op.m, op.cmul, op.mul)

*(op::LinearMap, b) = A_mul_B(op, b)

function A_mul_B{R,S}(op::LinearMap{R}, b::AbstractVector{S})
    A_mul_B!(Array(promote_type(R,S), op.m), op, b)
end
function A_mul_B{R,S}(op::LinearMap{R}, b::AbstractMatrix{S})
    A_mul_B!(Array(promote_type(R,S), op.m, size(b,2)), op, b)
end

A_mul_B!{T}(output, op::LinearMap{T,false}, b) = error("A*b not defined")
function A_mul_B!{T}(output, op::LinearMap{T,true}, b::AbstractVector)
    op.mul(output, b)
end
function A_mul_B!{T}(output, op::LinearMap{T,true}, b::AbstractMatrix)
    columns = [op.mul(output, b[:,i]) for i in 1:size(b,2)]
    hcat(columns...)
end

function Ac_mul_B{R,S}(op::LinearMap{R}, b::AbstractVector{S})
    Ac_mul_B!(Array(promote_type(R,S), op.n), op, b)
end
function Ac_mul_B{R,S}(op::LinearMap{R}, b::AbstractMatrix{S})
    Ac_mul_B!(Array(promote_type(R,S), op.n, size(b,2)), op, b)
end

Ac_mul_B!{T,A}(output, op::LinearMap{T,A,false}, b) = error("A'*b not defined")
function Ac_mul_B!{T,A}(output, op::LinearMap{T,A,true}, b::AbstractVector)
    op.cmul(output, b)
end
function Ac_mul_B!{T,A}(output, op::LinearMap{T,A,true}, b::AbstractMatrix)
    columns = [op.cmul(output, b[:,i]) for i in 1:size(b,2)]
    hcat(columns...)
end
