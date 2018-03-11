# function for arnoldi method

function arnoldi_procedure(A,m::Int64)
    # Procedure to generate orthogonal basis of the Krylov subspace
    n = size(A,1)
    Q = Matrix{Float64}(n,m+1)
    h = zeros(m+1,m)
    v = rand(n)
    w = Vector{Float64}(n)
    b=view(Q,:,1)
    copy!(b,v)
    normalize!(b)
    # fview(q,1)=normalize(V)
    for j = 1:m
        A_mul_B!(w,A,view(Q,:,j))

        # Gram-Schmidt Orthogonalization
        for i = 1:j
            q1 = view(Q,:,i)
            # h1 = view(h,i,j)
            h[i,j] = dot(w,q1)
            # zero dimensional subarray
            w = w-q1*h[i,j]
        end
        #h1 = view(h,j+1,j)
        #copy!(h1,norm(w,1))
        h[j+1,j] = norm(w,1)
        if h[j+1,j] == 0 # Found an Invariant Subspace
            return Q,h
        end
        q1=view(Q,:,j+1)
        copy!(q1,w/h[j+1,j])
    end
    return Q,h
end
