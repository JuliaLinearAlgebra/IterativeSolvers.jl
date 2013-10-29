import Base.LinAlg.BlasFloat

export eigvals_lanczos, svdvals_gkl

function lanczos{T<:BlasFloat}(K::KrylovSubspace{T})
    m = K.n
    αs = Array(T, m)
    βs = Array(T, m-1)
    for j=1:m-1
        w = nextvec(K)
        if j>1 w -= βs[j-1]*K.v[1] end
        w, y = orthogonalize(w, K, 1)
        αs[j] = y[1]
        βs[j] = norm(w)
        append!(K, w/βs[j])
    end
    αs[m]= dot(nextvec(K), lastvec(K))
    αs, βs
end

function eigvals_lanczos{T<:BlasFloat}(A::AbstractMatrix{T}, neigs::Int, verbose::Bool, tol::T, maxiter::Int)
    K = KrylovSubspace(A, 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    e1 = eigvals(SymTridiagonal(lanczos(K)...), 1, neigs)
    for iter=1:maxiter
        e0, e1 = e1, eigvals(SymTridiagonal(lanczos(K)...), 1, neigs)
        de = norm(e1-e0)
        if verbose println("Iteration ", iter, ": ", de) end
        if de < tol return e1 end
    end
    warn(string("Not converged: change in eigenvalues ", de, " exceeds specified tolerance of ", tol))
    e1
end
eigvals_lanczos{T<:BlasFloat}(A::AbstractMatrix{T}, neigs::Int=size(A,1), verbose=false, maxiter=size(A,1)) =eigvals_lanczos(A, neigs, verbose, size(A,1)^3*eps(T), maxiter)

#Golub-Kahan-Lanczos algorithm for the singular values
function svdvals_gkl{T<:BlasFloat}(A::AbstractMatrix{T})
    n = size(A, 2)
    α, β = Array(T, n), Array(T, n-1)
    K = KrylovSubspace(A, 1)
    initrand!(K)
    u, v = nextvec(K), lastvec(K)
    for k=1:n
        α[k]=norm(u)
        if k==n break end
        u/=α[k]
        v=A'*u-α[k]*v
        β[k]=norm(v)
        v/=β[k]
        u=A*v-β[k]*u
    end
    svdvals(Bidiagonal(α, β, true))
end

