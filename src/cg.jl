export cg, cg!

cg(A, b, Pl=1; kwargs...) =  cg!(zerox(A,b), A, b, Pl; kwargs...)

function cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(x, K, b, Pl; tol=tol, maxiter=maxiter)
end

function cg!(x, K::KrylovSubspace, b, Pl=1;
        tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2))
    resnorms = zeros(maxiter)

    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = isa(Pl, Function) ? Pl(r) : Pl\r
    γ = dot(r, z)
    for iter=1:maxiter
        append!(K, p)
        q = nextvec(K)
        α = γ/dot(p, q)
        # α>=0 || throw(PosSemidefException("α=$α"))
        update!(x, α, p)
        r -= α*q
        resnorms[iter] = norm(r)
        if resnorms[iter] < tol #Converged?
            resnorms = resnorms[1:iter]
            break
        end
        z = isa(Pl, Function) ? Pl(r) : Pl\r
        oldγ = γ
        γ = dot(r, z)
        β = γ/oldγ
        p = z + β*p
      end
    x, ConvergenceHistory(0<resnorms[end]<tol, tol, K.mvps, resnorms)
end

