export DGKS, ClassicalGramSchmidt, ModifiedGramSchmidt
export orthogonalize_and_normalize!

abstract OrthogonalizationMethod
immutable DGKS <: OrthogonalizationMethod end
immutable ClassicalGramSchmidt <: OrthogonalizationMethod end
immutable ModifiedGramSchmidt <: OrthogonalizationMethod end

# Default to MGS, good enough for solving linear systems.
@inline orthogonalize_and_normalize!{T}(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}) = orthogonalize_and_normalize!(V, w, h, ModifiedGramSchmidt)

function orthogonalize_and_normalize!{T}(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}, ::Type{DGKS})
    # Orthogonalize using BLAS-2 ops
    Ac_mul_B!(h, V, w)
    BLAS.gemv!('N', -one(T), V, h, one(T), w)
    nrm = norm(w)

    # Constant used by ARPACK.
    η = one(real(T)) / √2

    projection_size = norm(h)

    # Repeat as long as the DGKS condition is satisfied
    # Typically this condition is true only once.
    while nrm < η * projection_size
        correction = Ac_mul_B(V, w)
        projection_size = norm(correction)
        # w = w - V * correction
        BLAS.gemv!('N', -one(T), V, correction, one(T), w)
        @blas! h += one(T) * correction
        nrm = norm(w)
    end

    # Normalize; note that we already have norm(w).
    scale!(w, one(T) / nrm)

    nrm
end

function orthogonalize_and_normalize!{T}(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}, ::Type{ClassicalGramSchmidt})
    # Orthogonalize using BLAS-2 ops
    Ac_mul_B!(h, V, w)
    BLAS.gemv!('N', -one(T), V, h, one(T), w)
    nrm = norm(w)

    # Normalize
    scale!(w, one(T) / nrm)

    nrm
end

function orthogonalize_and_normalize!{T}(V::StridedVector{Vector{T}}, w::StridedVector{T}, h::StridedVector{T}, ::Type{ModifiedGramSchmidt})
    # Orthogonalize using BLAS-1 ops
    for i = 1 : length(V)
        h[i] = dot(V[i], w)
        @blas! w -= h[i] * V[i]
    end

    # Normalize
    nrm = norm(w)
    scale!(w, one(T) / nrm)

    nrm
end

function orthogonalize_and_normalize!{T}(V::StridedMatrix{T}, w::StridedVector{T}, h::StridedVector{T}, ::Type{ModifiedGramSchmidt})
    # Orthogonalize using BLAS-1 ops and column views.
    for i = 1 : size(V, 2)
        column = view(V, :, i)
        h[i] = dot(column, w)
        @blas! w -= h[i] * column
    end

    nrm = norm(w)
    scale!(w, one(T) / nrm)

    nrm
end