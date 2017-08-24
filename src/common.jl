import Base: A_ldiv_B!, \

using LinearMaps

export Identity

#### Type-handling
"""
    Adivtype(A, b)
Determine type of the division of an element of `b` against an element of `A`:
`typeof(one(eltype(b))/one(eltype(A)))`
"""
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))

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
solve!(out::AbstractArray{T},A::Int,b::AbstractArray{T}) where {T} = scale!(out,b, 1/A)
solve!(out::AbstractArray{T},A,b::AbstractArray{T}) where {T} = A_ldiv_B!(out,A,b)
solve!(out::AbstractArray{T},A::Function,b::AbstractArray{T}) where {T} = copy!(out,A(b))

# Identity preconditioner
struct Identity end

\(::Identity, x) = copy(x)
A_ldiv_B!(::Identity, x) = x
A_ldiv_B!(y, ::Identity, x) = copy!(y, x)