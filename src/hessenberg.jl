import Base.LinAlg: Givens, givensAlgorithm
import Base.A_ldiv_B!

struct FastHessenberg{T<:AbstractMatrix}
    H::T # H is assumed to be Hessenberg of size (m + 1) × m
end

@inline Base.size(H::FastHessenberg, args...) = size(H.H, args...)

"""
Solve Hy = rhs for a non-square Hessenberg matrix.
Note that `H` is also modified as is it converted
to an upper triangular matrix via Given's rotations
"""
function A_ldiv_B!(H::FastHessenberg, rhs)
    # Implicitly computes H = QR via Given's rotations
    # and then computes the least-squares solution y to
    # |Hy - rhs| = |QRy - rhs| = |Ry - Q'rhs|

    width = size(H, 2)

    # Hessenberg -> UpperTriangular; also apply to r.h.s.
    @inbounds for i = 1 : width
        c, s, _ = givensAlgorithm(H.H[i, i], H.H[i + 1, i])
        
        # Skip the first sub-diagonal since it'll be zero by design.
        H.H[i, i] = c * H.H[i, i] + s * H.H[i + 1, i]

        # Remaining columns
        @inbounds for j = i + 1 : width
            tmp = -conj(s) * H.H[i, j] + c * H.H[i + 1, j]
            H.H[i, j] = c * H.H[i, j] + s * H.H[i + 1, j]
            H.H[i + 1, j] = tmp
        end

        # Right hand side
        tmp = -conj(s) * rhs[i] + c * rhs[i + 1]
        rhs[i] = c * rhs[i] + s * rhs[i + 1]
        rhs[i + 1] = tmp
    end

    # Solve the upper triangular problem.
    U = UpperTriangular(view(H.H, 1 : width, 1 : width))
    A_ldiv_B!(U, view(rhs, 1 : width))
    nothing
end