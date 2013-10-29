import Base.LinAlg.BlasFloat

export eigvals_lanczos

function lanczos{T<:BlasFloat}(K::KrylovSubspace{T})
    m = K.n
    αs = Array(T, m)
    βs = Array(T, m-1)
    for j=1:m-1
        w = nextvec(K)
        αs[j] = dot(w, lastvec(K))
        orthogonalize!(w, K)
        βs[j] = norm(w)
        append!(K, w/βs[j])
    end
    αs[m]= dot(nextvec(K), lastvec(K))
    αs, βs
end

function eigvals_lanczos{T<:BlasFloat}(A::AbstractMatrix{T}, neigs::Int, verbose::Bool, tol::T, maxiter::Int)
    K = KrylovSubspace(A, 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    e0 = eigvals(SymTridiagonal(lanczos(K)...), 1, neigs)
    e1 = eigvals(SymTridiagonal(lanczos(K)...), 1, neigs)
    for iter=1:maxiter
        de = norm(e1-e0)
        if verbose println("Iteration ", iter, ": ", de) end
        if de < tol return e1 end
        e0, e1 = e1, eigvals(SymTridiagonal(lanczos(K)...), 1, neigs)
    end
    warn("Not converged")
    e1
end
eigvals_lanczos{T<:BlasFloat}(A::AbstractMatrix{T}, neigs::Int=size(A,1), verbose=false, maxiter=size(A,1)) =eigvals_lanczos(A, neigs, verbose, size(A,1)*size(A,2)*eps(T), maxiter)

