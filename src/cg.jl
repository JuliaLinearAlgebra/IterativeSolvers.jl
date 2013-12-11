export cg

function cg(A, b, x=nothing, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    K = KrylovSubspace(A, 1)
    if x==nothing || isempty(x)
        initrand!(K)
        x = lastvec(K)
    else
        init!(K, x)
    end
    cg(K, b, x, Pl; tol=tol, maxiter=maxiter)
end

function cg(K::KrylovSubspace, b, x, Pl=1; tol::Real=size(A,2)*eps(),
        maxIter::Integer=size(A,2))
    precondition(x) = isa(Pl, Function) ? Pl(x) : Pl\x
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

