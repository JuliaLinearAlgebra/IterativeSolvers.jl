export cg

function cg(A, b, Pl=x->x, x=nothing; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    K = KrylovSubspace(A, length(b), 1, Vector{eltype(b)}[])
    if x==nothing || isempty(x)
        initrand!(K)
        x = lastvec(K)
    else
        init!(K, x)
    end
    cg(K, b, Pl; tol=tol, maxiter=maxiter)
end

function cg(K::KrylovSubspace, b, Pl=x->x, x=nothing;
        tol::Real=size(A,2)*eps(), maxiter::Integer=size(A,2))
    precondition(x) = isa(Pl, Function) ? Pl(x) : Pl\x
    resnorms = zeros(maxiter)
 
    x = lastvec(K)
    r = b - nextvec(K)
    p = z = precondition(r)
    γ = dot(r, z)
    for iter=1:maxiter
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
    x, ConvergenceHistory(0<resnorms[end]<tol, tol, K.mvps, resnorms)
end

