using Compat

#----------------- Get the A matrix
function getDivGrad(n1,n2,n3)
    # the Divergence
    D1 = kron(speye(n3),kron(speye(n2),ddx(n1)))
    D2 = kron(speye(n3),kron(ddx(n2),speye(n1)))
    D3 = kron(ddx(n3),kron(speye(n2),speye(n1)))
    # DIV from faces to cell-centers

    Div = [D1 D2 D3]

    return Div*Div'
end
#----------------- 1D finite difference on staggered grid

# generate 1D derivatives
ddx(n)=spdiags(ones(n)*[-1 1],[0,1],n,n+1)

#------------- Build a diagonal matrix
function spdiags(B,d,m,n)
    # spdiags(B,d,m,n)
    # creates a sparse matrix from its diagonals
    d = d[:]
    p = length(d)

    len = zeros(p+1,1)
    for k = 1:p
        @compat len[k+1] = Int(len[k]+length(max(1,1-d[k]):min(m,n-d[k])))
    end
    @compat a = zeros(round(Int, len[p+1]), 3)
    for k = 1:p
        # Append new d[k]-th diagonal to compact form
        i = max(1,1-d[k]):min(m,n-d[k])
        @compat a[(round(Int, len[k])+1):round(Int, len[k+1]),:] = [i i+d[k] B[i+(m >= n)*d[k], k]]
    end

    @compat sparse(round(Int, a[:,1]), round(Int, a[:,2]), a[:,3], m, n)
end
