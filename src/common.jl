import  Base: eltype, length, ndims, real, size, *,
        A_mul_B!, Ac_mul_B, Ac_mul_B!

using   LinearMaps, Compat

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
    solve(A,b)
Solve `A\b` with a direct solver. When `A` is a function `A(b)` is dispatched instead.
"""
solve(A::Function,b) = A(b)

solve(A,b) = A\b

solve!{T}(out::AbstractArray{T},A::Int,b::AbstractArray{T}) = scale!(out,b, 1/A)

solve!{T}(out::AbstractArray{T},A,b::AbstractArray{T}) = A_ldiv_B!(out,A,b)
solve!{T}(out::AbstractArray{T},A::Function,b::AbstractArray{T}) = copy!(out,A(b))

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

# Non-allocating unsafe vector views
struct UnsafeVectorView{T} <: DenseVector{T}
    len::Int
    ptr::Ptr{T}
end

@inline UnsafeVectorView(parent::Vector, range::UnitRange) = UnsafeVectorView(length(range), pointer(parent) + sizeof(eltype(parent)) * (start(range) - 1))
@inline UnsafeVectorView(parent::Matrix, range::UnitRange, column::Int) = UnsafeVectorView(length(range), pointer(parent) + sizeof(eltype(parent)) * ((column - 1) * size(parent, 1) + start(range) - 1))
@inline Base.size(v::UnsafeVectorView) = (v.len,)
@inline Base.getindex(v::UnsafeVectorView, idx) = unsafe_load(v.ptr, idx)
@inline Base.setindex!(v::UnsafeVectorView, value, idx) = unsafe_store!(v.ptr, value, idx)
@inline Base.length(v::UnsafeVectorView) = v.len
@inline Base.pointer(v::UnsafeVectorView{T}) where {T} = v.ptr
@compat Base.IndexStyle(::Type{V}) where {V <: UnsafeVectorView} = IndexLinear()

# Identity preconditioner
immutable Identity end

Base.:\(::Identity, x) = copy(x)
Base.A_ldiv_B!(::Identity, x) = x
Base.A_ldiv_B!(y, ::Identity, x) = copy!(y, x)