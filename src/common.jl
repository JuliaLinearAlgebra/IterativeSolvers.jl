import LinearAlgebra: ldiv!, \

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
zerox(A, b::GPUArray) = zerox(A::GPUArray, b) = throw("Please pre-allocate the result vector on the GPU and use an inplace function.")

"""
No-op preconditioner
"""
struct Identity end

\(::Identity, x) = copy(x)
ldiv!(::Identity, x) = x
ldiv!(y, ::Identity, x) = copyto!(y, x)
