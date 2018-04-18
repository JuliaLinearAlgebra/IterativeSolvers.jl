import LinearAlgebra: A_ldiv_B!, \

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
zerox(A, b) = zeros(Adivtype(A, b), size(A, 2))

"""
No-op preconditioner
"""
struct Identity end

\(::Identity, x) = copy(x)
A_ldiv_B!(::Identity, x) = x
A_ldiv_B!(y, ::Identity, x) = copyto!(y, x)
