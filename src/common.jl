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

Build a zeros vector similar to `b` of eltype `Adivtype(A,b)`.
"""
function zerox(A, b)
    T = Adivtype(A, b)
    x = similar(b, T, size(A, 2))
    fill!(x, zero(T))
    return x
end

"""
No-op preconditioner
"""
struct Identity end

\(::Identity, x) = copy(x)
ldiv!(::Identity, x) = x
ldiv!(y, ::Identity, x) = copyto!(y, x)

"""
Conjugated and unconjugated dot products
"""
# Conjugated dot product
_dot(x, ::Val{true}) = sum(abs2, x)  # for x::Complex, returns Real
_dot(x, y, ::Val{true}) = dot(x, y)

# Unconjugated dot product
_dot(x, ::Val{false}) = sum(xₖ^2 for xₖ in x)
_dot(x, y, ::Val{false}) = sum(prod, zip(x,y))
