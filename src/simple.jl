#Simple methods

export ev_power, ev_ii

#Power method for finding largest eigenvalue and its eigenvector
function ev_power{T<:Number}(K::KrylovSubspace{T}, maxiter::Int, tol::T)
    θ = de = 0
    v = Array(T, K.n)
    for k=1:maxiter
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        de = norm(y - θ*v)
        if de <= tol * abs(θ) return θ, v end
        appendunit!(K, y)
    end
    warn(string("Not converged: change in eigenvector ", de, " exceeds specified tolerance of ", tol))
    θ, v
end

function ev_power{T<:Number}(A::AbstractMatrix{T}, maxiter::Int, tol::T)
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_power(K, maxiter, tol)
end
ev_power{T<:BlasFloat}(A::AbstractMatrix{T}, maxiter::Int=size(A,1)) = ev_power(A, maxiter, eps(T))

#Inverse iteration/inverse power method
function ev_ii{T<:Number}(K::KrylovSubspace{T}, σ::T, maxiter::Int, tol::T)
    θ = de = 0
    v = Array(T, K.n)
    for k=1:maxiter
        v = lastvec(K)
        y = (K.A-σ*diagm(ones(K.n)))\v
        θ = dot(v, y)
        de = norm(y - θ*v)
        if de <= tol * abs(θ) return σ+1/θ, y/θ end
        appendunit!(K, y)
    end
    warn(string("Not converged: change in eigenvector ", de, " exceeds specified tolerance of ", tol))
    σ+1/θ, y/θ
end

function ev_ii{T<:Number}(A::AbstractMatrix{T}, σ::T, maxiter::Int, tol::T)
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_ii(K, σ, maxiter, tol)
end
ev_ii{T<:BlasFloat}(A::AbstractMatrix{T}, σ::T, maxiter::Int=size(A,1)) = ev_ii(A, σ, maxiter, eps(T))

