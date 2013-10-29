#Simple methods

export ev_power, ev_ii, ev_rqi

#Power method for finding largest eigenvalue and its eigenvector
function ev_power{T<:Number}(K::KrylovSubspace{T}, maxiter::Int, tol::T)
    θ = de = 0
    v = Array(T, K.n)
    for k=1:maxiter
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        de = norm(y - θ*v)
        if de <= tol * abs(θ) return Eigenpair(θ, v) end
        appendunit!(K, y)
    end
    warn(string("Not converged: change in eigenvector ", de, " exceeds specified tolerance of ", tol))
    Eigenpair(θ, v)
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
        y = (K.A-σ*eye(K))\v
        θ = dot(v, y)
        de = norm(y - θ*v)
        if de <= tol * abs(θ) return Eigenpair(σ+1/θ, y/θ) end
        appendunit!(K, y)
    end
    warn(string("Not converged: change in eigenvector ", de, " exceeds specified tolerance of ", tol))
    Eigenpair(σ+1/θ, y/θ)
end

function ev_ii{T<:Number}(A::AbstractMatrix{T}, σ::T, maxiter::Int, tol::T)
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_ii(K, σ, maxiter, tol)
end
ev_ii{T<:BlasFloat}(A::AbstractMatrix{T}, σ::T, maxiter::Int=size(A,1)) = ev_ii(A, σ, maxiter, eps(T))

#Rayleigh quotient iteration
#XXX Doesn't work well
function ev_rqi{T<:Number}(K::KrylovSubspace{T}, σ::T, maxiter::Int, tol::T)
    v = lastvec(K)
    ρ = dot(v, nextvec(K))
    for k=1:maxiter
        y = (K.A-ρ*eye(K))\v
        θ = norm(y)
        ρ += dot(y,v)/θ^2
        v = y/θ
        if θ >= 1/tol break end
    end
    Eigenpair(ρ, v)
end
function ev_rqi{T<:Number}(A::AbstractMatrix{T}, σ::T, maxiter::Int, tol::T)
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_rqi(K, σ, maxiter, tol)
end
ev_rqi{T<:BlasFloat}(A::AbstractMatrix{T}, σ::T, maxiter::Int=size(A,1)) = ev_rqi(A, σ, maxiter, eps(T))


