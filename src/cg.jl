export cg

function cg(A, b, x=nothing; tol=size(A,2)*eps(), maxIter=size(A,2),
        Preconditioner=1)
    K = KrylovSubspace(A, 1)
    if x==nothing || isempty(x)
        initrand!(K)
        x = lastvec(K)
    else
        init!(K, x)
    end
    cg(K, b, x; tol=tol, maxIter=maxIter, Preconditioner=Preconditioner)
end

function cg(K::KrylovSubspace, b, x; tol::Real=size(A,2)*eps(),
        maxIter::Integer=size(A,2), Preconditioner=1)
    precondition(x) = isa(Preconditioner, Function) ? Preconditioner(x) :
        Preconditioner\x
    resnorms = zeros(maxIter)
 
    r = b - nextvec(K)
    p = z = precondition(r)
    γ = dot(r, z)
    for iter=1:maxIter
        append!(K, p)
        q = nextvec(K)
        α = γ/dot(p, q)
        α>=0 || throw(PosSemidefException("α=$α"))
        x += α*p
        r -= α*q

        resnorms[iter] = norm(r)
        if resnorms[iter] < tol #Converged?
            resnorms = resnorms[1:iter]
            break
        end
        z = precondition(r)
        oldγ = γ
        γ = dot(r, z)
        β = γ/oldγ
        p = z + β*p
      end
    x, ConvergenceHistory(0<resnorms[end]<tol, tol, resnorms)
end

