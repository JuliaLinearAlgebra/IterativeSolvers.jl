export cg, cg!

cg(A, b, Pl=1; kwargs...) =  cg!(zerox(A,b), A, b, Pl; kwargs...)

function cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)

    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg_method!(history, x, K, b, Pl; tol=tol, maxiter=maxiter)
    x, history
end

function cg_method!(log::ConvergenceHistory, x, K::KrylovSubspace, b, Pl=1;
        tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2))

    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = Pl\r
    γ = dot(r, z)
    for iter=1:maxiter
        nextiter!(log)
        append!(K, p)
        q = nextvec(K)
        α = γ/dot(p, q)
        # α>=0 || throw(PosSemidefException("α=$α"))
        update!(x, α, p)
        r -= α*q
        resnorm = norm(r)
        push!(log,:resnorm,resnorm)
        resnorm < tol && break
        z = Pl\r
        oldγ = γ
        γ = dot(r, z)
        β = γ/oldγ
        p = z + β*p
    end
    shrink!(log)
    setmvps(log, K.mvps)
    setconv(log, 0<=norm(r)<tol)
    x
end
