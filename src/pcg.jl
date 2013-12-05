
function pcg(A,b,tol=1e-2,maxIter=100,M=1,x=[])
    resvec = zeros(maxIter)
    Af(x) =  isa(A,Function) ? A(x) : A*x
    Mf(x) =  isa(M,Function) ? M(x) : M\x

    if isempty(x)
        x = zeros(size(b,1))
        r = b
    else
       r = b - Af(x)
    end
    z = Mf(r)
    p = z
    nr0  = norm(b)
    for i=1:maxIter
        Ap = Af(p)
        gamma = dot(r,z)
        alpha = gamma/dot(p,Ap)

        x += alpha*p
        r -= alpha*Ap
        resvec[i]  = norm(r)/nr0
        if resvec[i] < tol
            return x,0,resvec[i],i,resvec[1:i]
        end
        z = Mf(r)

        beta = dot(z,r)/gamma[1]
        p = z + beta*p
     end
    return x,-1,resvec[maxIter],maxIter,resvec
 end

export pcg
