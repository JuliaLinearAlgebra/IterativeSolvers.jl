export DGKS, ClassicalGramSchmidt, ModifiedGramSchmidt
export orthogonalize_and_normalize!

abstract type OrthogonalizationMethod end
struct DGKS <: OrthogonalizationMethod end
struct ClassicalGramSchmidt <: OrthogonalizationMethod end
struct ModifiedGramSchmidt <: OrthogonalizationMethod end

# Default to MGS, good enough for solving linear systems.
@inline orthogonalize_and_normalize!(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}) where {T} = 
	orthogonalize_and_normalize!(V, w, h, ModifiedGramSchmidt())

function orthogonalize_and_normalize!(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}, ::DGKS) where {T}
    # Orthogonalize using BLAS-2 ops
    mul!(h, adjoint(V), w)
    mul!(w, V, h, -one(T), one(T))
    nrm = norm(w)

    # Constant used by ARPACK.
    η = one(real(T)) / √2

    projection_size = norm(h)

    # Repeat as long as the DGKS condition is satisfied
    # Typically this condition is true only once.
    while nrm < η * projection_size
        correction = adjoint(V)*w
        projection_size = norm(correction)
        # w = w - V * correction
        mul!(w, V, correction, -one(T), one(T))
        h .+= correction
        nrm = norm(w)
    end

    # Normalize; note that we already have norm(w).
    w .*= inv(nrm)

    nrm
end

function orthogonalize_and_normalize!(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}, ::ClassicalGramSchmidt) where {T}
    # Orthogonalize using BLAS-2 ops
    mul!(h, adjoint(V), w)
    mul!(w, V, h, -one(T), one(T))
    nrm = norm(w)

    # Normalize
    w .*= inv(nrm)

    nrm
end

function orthogonalize_and_normalize!(V::StridedVector{Vector{T}}, w::StridedVector{T}, h::StridedVector{T}, ::ModifiedGramSchmidt) where {T}
    # Orthogonalize using BLAS-1 ops
    for i = 1 : length(V)
        h[i] = dot(V[i], w)
        w .-= h[i] .* V[i]
    end

    # Normalize
    nrm = norm(w)
    w .*= inv(nrm)

    nrm
end

function orthogonalize_and_normalize!(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}, ::ModifiedGramSchmidt) where {T}
    # Orthogonalize using BLAS-1 ops and column views.
    for i = 1 : size(V, 2)
        column = view(V, :, i)
        h[i] = dot(column, w)
        w .-= h[i] .* column
    end

    nrm = norm(w)
    w .*= inv(nrm)

    nrm
end
