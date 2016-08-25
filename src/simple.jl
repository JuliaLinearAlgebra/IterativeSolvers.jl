#Simple methods

export eigvals_power, eigvals_ii, eigvals_rqi


function powm_method!{T}(log::ConvergenceHistory, K::KrylovSubspace{T};
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n
    )
    θ = zero(T)
    v = Array(T, K.n)
    for iter=1:maxiter
        nextiter!(log,mvps=1)
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        resnorm = real(norm(y - θ*v))
        push!(log, :resnorm, resnorm)
        resnorm <= tol*abs(θ) && (setconv(log, resnorm >= 0); break)
        appendunit!(K, y)
    end
    Eigenpair(θ, v)
end

function eigvals_power(A, x=nothing; tol::Real=size(A,1)^3*eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    eigvals_power(K; tol=tol, maxiter=maxiter)
end
function eigvals_power(K::KrylovSubspace, x=nothing;
    tol::Real=size(A,1)^3*eps(), maxiter::Int=size(A,1)
    )
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    pair = powm_method!(history, K; tol=tol, maxiter=maxiter)
    pair, history
end

#Inverse iteration/inverse power method
function invpowm_method!{T}(log::ConvergenceHistory, K::KrylovSubspace{T}, σ::Number=zero(T);
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n
    )
    θ = zero(T)
    v = Array(T, K.n)
    y = Array(T, K.n)
    σ = convert(T, σ)
    for iter=1:maxiter
        nextiter!(log,mvps=1)
        v = lastvec(K)
        y = (K.A-σ*eye(K))\v
        θ = dot(v, y)
        resnorm = norm(y - θ*v)
        push!(log, :resnorm, resnorm)
        resnorm <= tol*abs(θ) && (setconv(log, resnorm >= 0); break)
        appendunit!(K, y)
    end
    shrink!(log)
    Eigenpair(σ+1/θ, y/θ)
end

function eigvals_ii(A, σ::Number, x=nothing; tol::Real=eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    eigvals_ii(K, σ; tol=tol, maxiter=maxiter)
end
function eigvals_ii(K::KrylovSubspace, σ::Number, x=nothing;
    tol::Real=size(A,1)^3*eps(), maxiter::Int=size(A,1)
    )
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    pair = invpowm_method!(history, K, σ; tol=tol, maxiter=maxiter)
    pair, history
end

#Rayleigh quotient iteration
#XXX Doesn't work well
function rqi_method!{T}(log::ConvergenceHistory, K::KrylovSubspace{T}, σ::Number; tol::Real=eps(T), maxiter::Int=K.n)
    v = lastvec(K)
    ρ = dot(v, nextvec(K))
    for iter=1:maxiter
        nextiter!(log,mvps=1)
        y = (K.A-ρ*eye(K))\v
        θ = norm(y)
        ρ += dot(y,v)/θ^2
        v = y/θ
        resnorm=1/θ
        push!(log,:resnorm,resnorm)
        θ >= 1/tol && (setconv(log, resnorm >= 0); break)
    end
    shrink!(history)
    Eigenpair(ρ, v)
end

function eigvals_rqi(A, σ::Number, x=nothing; tol::Real=eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    eigvals_rqi(K, σ; tol=tol, maxiter=maxiter)
end
function eigvals_rqi{T}(K::KrylovSubspace{T}, σ::Number, x=nothing; tol::Real=eps(T), maxiter::Int=K.n)
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    rqi_method(history, K, σ, tol=tol, maxiter=maxiter)
    x, history
end
