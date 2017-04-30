import  Base: eltype, eps, length, ndims, real, size, *, \,
        A_mul_B!, Ac_mul_B, Ac_mul_B!

import LinearMaps.FunctionMap

export  A_mul_B

using   LinearMaps

# Improve readability of iterative methods
\(f::Function, b) = f(b)
*(f::Function, b) = f(b)

#### Type-handling
"""
    Adivtype(A, b)
If A is a function map then:
`typeof(one(eltype(b)))`
Determine type of the division of an element of `b` against an element of `A`:
`typeof(one(eltype(b))/one(eltype(A)))`
"""
Adivtype(A, b) = isa(A, FunctionMap) ? typeof(one(eltype(b))) : typeof(one(eltype(b))/one(eltype(A)))

"""
    Amultype(A, x)
If A is a function map then:
`typeof(one(eltype(b)))`
Determine type of the multiplication of an element of `b` with an element of `A`:
`typeof(one(eltype(A))*one(eltype(x)))`
"""
Amultype(A, x) = isa(A, FunctionMap) ? typeof(one(eltype(b))) :typeof(one(eltype(A))*one(eltype(x)))

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
