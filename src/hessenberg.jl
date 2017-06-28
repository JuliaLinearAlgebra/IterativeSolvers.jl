import Base.LinAlg: Givens, givensAlgorithm

mutable struct MyHessenberg{T<:AbstractMatrix}
    H::T
end

@inline Base.size(H::MyHessenberg, args...) = size(H.H, args...)

function Base.A_mul_B!(G::Givens, H::MyHessenberg)
    m, n = size(H)
    @inbounds for i = G.i1 : n
        a1, a2 = H.H[G.i1, i], H.H[G.i2, i]
        H.H[G.i1, i] =       G.c  * a1 + G.s * a2
        H.H[G.i2, i] = -conj(G.s) * a1 + G.c * a2
    end
    return H
end

function solve!(H::MyHessenberg, rhs)
    width = size(H, 2)

    # MyHessenberg -> UpperTriangular; also apply to r.h.s.
    @inbounds for i = 1 : width
        c, s, r = givensAlgorithm(H.H[i, i], H.H[i + 1, i])
        
        # Diagonal element
        H.H[i, i] = c * H.H[i, i] + s * H.H[i + 1, i]

        # Remaining columns
        @inbounds for j = i + 1 : width
            tmp = -s * H.H[i, j] + c * H.H[i + 1, j]
            H.H[i, j] = c * H.H[i, j] + s * H.H[i + 1, j]
            H.H[i + 1, j] = tmp
        end

        # Right hand side
        tmp = -s * rhs[i] + c * rhs[i + 1]
        rhs[i] = c * rhs[i] + s * rhs[i + 1]
        rhs[i + 1] = tmp
    end

    # Solve the upper triangular problem.
    U = UpperTriangular(@view(H.H[1 : width, 1 : width]))
    U \ @view(rhs[1 : width])
end