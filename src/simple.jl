#Simple methods

export powm, invpowm, rqi

function powm(A; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    powm(K; kwargs...)
end

powm(K::KrylovSubspace; kwargs...) = powm_method(K; kwargs...)

function powm(::Type{Master}, A; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    powm(Master, K; kwargs...)
end

function powm{T}(::Type{Master}, K::KrylovSubspace{T};
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    eigs = powm_method(K; verbose=verbose,log=log,tol=tol,maxiter=maxiter)
    shrink!(log)
    plot && showplot(log)
    eigs, log
end

#Power method for finding largest eigenvalue and its eigenvector
function powm_method{T}(K::KrylovSubspace{T};
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    verbose::Bool=false, log::MethodLog=DummyHistory()
    )
    θ = zero(T)
    v = Array(T, K.n)
    for iter=1:maxiter
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        resnorm = real(norm(y - θ*v))
        nextiter!(log)
        push!(log, :resnorm, resnorm)
        resnorm <= tol*abs(θ) && (setconv(log, resnorm >= 0); break)
        appendunit!(K, y)
    end
    setmvps(log, K.mvps)
    Eigenpair(θ, v)
end

function invpowm(A; shift::Number=0, x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    invpowm(K; shift=shift, kwargs...)
end

invpowm(K::KrylovSubspace; shift::Number=0, kwargs...) = invpowm_method(K, shift; kwargs...)

function invpowm(::Type{Master}, A; shift::Number=0, x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    invpowm(Master, K; shift=shift, kwargs...)
end

function invpowm{T}(::Type{Master}, K::KrylovSubspace{T};
    shift::Number=0, tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    eigs = invpowm_method(K, shift;
        verbose=verbose, log=log, tol=tol, maxiter=maxiter
        )
    shrink!(log)
    plot && showplot(log)
    eigs, log
end

#Inverse iteration/inverse power method
function invpowm_method{T}(K::KrylovSubspace{T}, σ::Number;
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    log::MethodLog=DummyHistory(), verbose::Bool=false
    )
    θ = zero(T)
    v = Array(T, K.n)
    y = Array(T, K.n)
    σ = convert(T, σ)
    for iter=1:maxiter
        v = lastvec(K)
        y = (K.A-σ*eye(K))\v
        θ = dot(v, y)
        resnorm = norm(y - θ*v)
        nextiter!(log)
        push!(log, :resnorm, resnorm)
        resnorm <= tol*abs(θ) && (setconv(log, resnorm >= 0); break)
        appendunit!(K, y)
    end
    setmvps(log, K.mvps)
    Eigenpair(σ+1/θ, y/θ)
end

function rqi(A; shift::Number=0, x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    rqi(K; shift=shift, kwargs...)
end

rqi(K::KrylovSubspace; shift::Number=0, kwargs...) = rqi_method(K, shift; kwargs...)

function rqi(::Type{Master}, A; shift::Number=0, x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    rqi(Master, K; shift=shift, kwargs...)
end

function rqi{T}(::Type{Master}, K::KrylovSubspace{T};
    shift::Number=0, tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    eigs = rqi_method(K, shift; verbose=verbose, log=log, tol=tol, maxiter=maxiter)
    shrink!(log)
    plot && showplot(log)
    eigs, log
end

#Rayleigh quotient iteration
#XXX Doesn't work well
function rqi_method{T}(K::KrylovSubspace{T}, σ::Number;
    x=nothing, tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    log::MethodLog=DummyHistory(), verbose::Bool=false
    )
    v = lastvec(K)
    ρ = dot(v, nextvec(K))
    resnorms=zeros(eltype(real(one(T))), maxiter)
    for iter=1:maxiter
        y = (K.A-ρ*eye(K))\v
        θ = norm(y)
        ρ += dot(y,v)/θ^2
        v = y/θ
        resnorm=1/θ
        nextiter!(log)
        push!(log,:resnorm,resnorm)
        θ >= 1/tol && (setconv(log, resnorm >= 0); break)
    end
    setmvps(log, K.mvps)
    Eigenpair(ρ, v)
end
