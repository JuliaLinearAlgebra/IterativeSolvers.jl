import Base.LinAlg: Givens, givensAlgorithm

mutable struct Hessenberg{T<:AbstractMatrix}
    H::T # H is assumed to be Hessenberg of size (m + 1) × m
end

@inline Base.size(H::Hessenberg, args...) = size(H.H, args...)

function solve!(H::Hessenberg, rhs)
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
    U = UpperTriangular(view(H.H, 1 : width, 1 : width))
    U \ UnsafeVectorView(rhs, 1 : width)
end