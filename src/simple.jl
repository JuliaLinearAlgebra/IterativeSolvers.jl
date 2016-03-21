#Simple methods

export eigvals_power, eigvals_ii, eigvals_rqi

#Power method for finding largest eigenvalue and its eigenvector
function eigvals_power{T}(K::KrylovSubspace{T}; tol::Real=eps(T)*K.n^3, maxiter::Int=K.n)
    θ = zero(T)
    v = Array(T, K.n)
    resnorms=zeros(eltype(real(one(T))), maxiter)
    for iter=1:maxiter
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        resnorms[iter] = real(norm(y - θ*v))
        if resnorms[iter] <= tol * abs(θ)
            resnorms=resnorms[1:iter]
            break
        end
        appendunit!(K, y)
    end
    Eigenpair(θ, v), ConvergenceHistory(0<resnorms[end]<tol*abs(θ), tol, K.mvps, resnorms)

end

function eigvals_power(A, x=nothing; tol::Real=size(A,1)^3*eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    eigvals_power(K; tol=tol, maxiter=maxiter)
end

#Inverse iteration/inverse power method
function eigvals_ii{T}(K::KrylovSubspace{T}, σ::Number=zero(T); tol::Real=eps(T)*K.n^3, maxiter::Int=K.n)
    θ = zero(T)
    v = Array(T, K.n)
    y = Array(T, K.n)
    σ = convert(T, σ)
    resnorms=zeros(eltype(real(one(T))), maxiter)
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
    Eigenpair(σ+1/θ, y/θ), ConvergenceHistory(0<resnorms[end]<tol*abs(θ), tol, K.mvps, resnorms)
end

function eigvals_ii(A, σ::Number, x=nothing; tol::Real=eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    eigvals_ii(K, σ; tol=tol, maxiter=maxiter)
end

#Rayleigh quotient iteration
function eigvals_rqi{T}(K::KrylovSubspace{T}, σ::Number=zero(T); tol::Real=eps(T)*K.n^3, maxiter::Int=K.n)
    θ = zero(T)
    v = Array(T, K.n)
    y = Array(T, K.n)
    σ = convert(T, σ)
    resnorms=zeros(eltype(real(one(T))), maxiter)
    for iter=1:maxiter
        v = lastvec(K)
        y = (K.A-σ*eye(K))\v
        θ = dot(v, y)
	σ = σ+1/θ
        resnorms[iter] = norm(y - θ*v)
        if resnorms[iter] <= tol * abs(θ)
            resnorms=resnorms[1:iter]
            break
        end
        appendunit!(K, y)
    end
    Eigenpair(σ, y/θ), ConvergenceHistory(0<resnorms[end]<tol*abs(θ), tol, K.mvps, resnorms)
end
function eigvals_rqi(A, σ::Number, x=nothing; tol::Real=eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    eigvals_rqi(K, σ; tol=tol, maxiter=maxiter)
end
