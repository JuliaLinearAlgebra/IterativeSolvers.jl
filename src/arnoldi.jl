# function for arnoldi method 

function arnoldi_procedure(A::StridedMatrix{T}, w::StridedVector{T}, h::StridedMatrix{T}, v::StridedVector{T},m::Int64) where {T}
    # Procedure to generate orthogonal basis of the Krylov subspace
    q = Matrix(m,size(v))
    q[1,:] = normalize(v)
    for j = 1:m
        A_mul_B!(w,A,q[j,:])
        
        # Gram-Schmidt Orthogonalization
        for i = 1:j
            h[i,j] = A_mul_B(w',q[i,:])
            w = w-q[i,:]*h[i,j]     
        end
        h[j+1,j] = norm(w,1)
        if h[j+1,j] == 0 # Found an Invariant Subspace
            return q,h
        end
        q[j+1,:] = w/h[j+1,j]
    end
    return q,h
end
