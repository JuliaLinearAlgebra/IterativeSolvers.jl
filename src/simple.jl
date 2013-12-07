#Simple methods

export ev_power, ev_ii, ev_rqi

#Power method for finding largest eigenvalue and its eigenvector
function ev_power{T}(K::KrylovSubspace{T}, maxiter::Int=K.n, tol::Real=eps(T)*K.n^3)
    θ = zero(T)
    v = Array(T, K.n)
    resnorms=zeros(maxiter)
    for iter=1:maxiter
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        resnorms[iter] = norm(y - θ*v)
        if resnorms[iter] <= tol * abs(θ)
            resnorms=resnorms[1:iter]
            break
        end
        appendunit!(K, y)
    end
    Eigenpair(θ, v), ConvergenceHistory(0<resnorms[end]<tol, tol, resnorms)

end

function ev_power(A, maxiter::Int=size(A,1), tol::Real=size(A,1)^3*eps())
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_power(K, maxiter, tol)
end

#Inverse iteration/inverse power method
function ev_ii{T}(K::KrylovSubspace{T}, σ::Number=zero(T), maxiter::Int=K.n, tol::Real=eps(T)*K.n^3)
    θ = zero(T)
    v = Array(T, K.n)
    y = Array(T, K.n)
    resnorms=zeros(maxiter)
    for iter=1:maxiter
        v = lastvec(K)
        y = (K.A-σ*eye(K))\v
        θ = dot(v, y)
        resnorms[iter] = norm(y - θ*v)
        if resnorms[iter] <= tol * abs(θ)
            resnorms=resnorms[1:iter]
            break
        end
        appendunit!(K, y)
    end
    Eigenpair(σ+1/θ, y/θ), ConvergenceHistory(0<resnorms[end]<tol, tol, resnorms)
end

function ev_ii(A, σ::Number, maxiter::Int=size(A,1), tol::Real=eps())
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_ii(K, σ, maxiter, tol)
end

#Rayleigh quotient iteration
#XXX Doesn't work well
function ev_rqi(K::KrylovSubspace, σ::Number, maxiter::Int, tol::Real)
    v = lastvec(K)
    ρ = dot(v, nextvec(K))
    resnorms=zeros(maxiter)
    for iter=1:maxiter
        y = (K.A-ρ*eye(K))\v
        θ = norm(y)
        ρ += dot(y,v)/θ^2
        v = y/θ
        resnorms[iter]=1/θ
        if θ >= 1/tol 
            resnorms=resnorms[1:iter]
            break
        end
    end
    Eigenpair(ρ, v), ConvergenceHistory(0<resnorms[end]<tol, tol, resnorms)
end
function ev_rqi(A, σ::Number, maxiter::Int=size(A,1), tol::Real=eps())
    K = KrylovSubspace(A, 1)
    initrand!(K)
    ev_rqi(K, σ, maxiter, tol)
end

