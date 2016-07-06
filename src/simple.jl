#Simple methods

export powm, master_powm, invpowm, master_invpowm, rqi, master_rqi

function powm(A; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    powm(K; kwargs...)
end

function powm(K::KrylovSubspace; kwargs...)
    eigs, _ = powm_method(K; kwargs...)
    eigs
end

function powm(::Type{Master}, A; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    powm(Master, K; kwargs...)
end

function powm{T}(::Type{Master}, K::KrylovSubspace{T};
    tol::Real=epsMaster, (T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false
    )
    log = MethodLog(maxiter)
    add!(log,:resnorm)
    eigs, conv = powm_method(K; verbose=verbose,log=log,tol=tol,maxiter=maxiter)
    shrink!(log)
    plot && showplot(log)
    eigs, ConvergenceHistory(conv, tol, K.mvps, log)
end

#Power method for finding largest eigenvalue and its eigenvector
function powm_method{T}(K::KrylovSubspace{T};
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    verbose::Bool=false, log::MethodLog=MethodLog()
    )
    θ = zero(T)
    v = Array(T, K.n)
    for iter=1:maxiter
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        resnorm = real(norm(y - θ*v))
        next!(log)
        push!(log, :resnorm, resnorm)
        resnorm <= tol*abs(θ) && break
        appendunit!(K, y)
    end
    Eigenpair(θ, v), isconverged(log,:resnorm,tol*abs(θ))
end

function invpowm(A, σ::Number; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    invpowm(K, σ; kwargs...)
end

function invpowm(K::KrylovSubspace, σ::Number; kwargs...)
    eigs, _ = invpowm_method(K, σ; kwargs...)
    eigs
end

function invpowm(::Type{Master}, A, σ::Number; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    invpowm(Master, K, σ; kwargs...)
end

function invpowm{T}(::Type{Master}, K::KrylovSubspace{T}, σ::Number;
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false
    )
    log = MethodLog(maxiter)
    add!(log,:resnorm)
    eigs, conv = invpowm_method(K, σ;
        verbose=verbose, log=log, tol=tol, maxiter=maxiter
        )
    shrink!(log)
    plot && showplot(log)
    eigs, ConvergenceHistory(conv, tol, K.mvps, log)
end

#Inverse iteration/inverse power method
function invpowm{T}(::Type{Master}, K::KrylovSubspace{T}, σ::Number;
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    log::MethodLog=MethodLog(), verbose::Bool=false
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
        next!(log)
        push!(log, :resnorm, resnorm)
        resnorm <= tol*abs(θ) && break
        appendunit!(K, y)
    end
    Eigenpair(σ+1/θ, y/θ), isconverged(log,:resnorm,tol*abs(θ))
end

function rqi(A, σ::Number; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    rqi(K, σ; kwargs...)
end

function rqi(K::KrylovSubspace, σ::Number; kwargs...)
    eigs, _ = rqi_method(K, σ; kwargs...)
    eigs
end

function rqi(::Type{Master}, A, σ::Number; x=nothing, kwargs...)
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))
    rqi(Master, K, σ; kwargs...)
end

function rqi{T}(::Type{Master}, K::KrylovSubspace{T}, σ::Number;
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false
    )
    log = MethodLog(maxiter)
    add!(log,:resnorm)
    eigs, conv = rqi_method(K, σ; verbose=verbose, log=log, tol=tol, maxiter=maxiter)
    shrink!(log)
    plot && showplot(log)
    eigs, ConvergenceHistory(conv, tol, K.mvps, log)
end

#Rayleigh quotient iteration
#XXX Doesn't work well
function rqi_method{T}(K::KrylovSubspace{T}, σ::Number;
    x=nothing, tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    log::MethodLog=MethodLog(), verbose::Bool=false
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
        next!(log)
        push!(log,:resnorm,resnorm)
        θ >= 1/tol && break
    end
    Eigenpair(ρ, v), isconverged(log,:resnorm,tol)
end
