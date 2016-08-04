#Simple methods
export powm, invpowm, rqi

####################
# API method calls #
####################

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
    eigs..., log
end

#########################
# Method Implementation #
#########################

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
    θ, v
end

####################
# API method calls #
####################

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
    eigs..., log
end

#########################
# Method Implementation #
#########################

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
    σ+1/θ, y/θ
end

####################
# API method calls #
####################

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
    eigs..., log
end

#########################
# Method Implementation #
#########################

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
    ρ, v
end

#################
# Documentation #
#################

#Initialize parameters
doc1_call = """    powm(A)
    powm(K)
    powm(Master, A)
    powm(Master, K)
"""
doc2_call = """    invpowm(A)
    invpowm(K)
    invpowm(Master, A)
    invpowm(Master, K)
"""
doc3_call = """    rqi(A)
    rqi(K)
    rqi(Master, A)
    rqi(Master, K)
"""
doc1_msg = """Find biggest eigenvalue of `A` and its associated eigenvector
using the power method.
"""
doc2_msg = """Find closest eigenvalue of `A` to `shift` and its associated eigenvector
using the inverse power iteration method.
"""
doc3_msg = """Try find closest eigenvalue of `A` to `shift` and its associated eigenvector
using the rayleigh quotient iteration method. This method converges rapidly
but is not guaranteed to compute the eigenvalue closes to `shift`.
"""
doc1_karg = ""
doc2_karg = "* `shift::Number=0`: shift to be applied to matrix A."
doc3_karg = "* `shift::Number=0`: shift to be applied to matrix A."

doc1_version = (powm, doc1_call, doc1_msg, doc1_karg)
doc2_version = (invpowm, doc2_call, doc2_msg, doc2_karg)
doc3_version = (rqi, doc3_call, doc3_msg, doc3_karg)

#Build docs
for (func, call, msg, karg) in [doc1_version, doc2_version, doc3_version]
@doc """
$call

$msg

If [`Master`](@ref) is given, method will output a tuple `x, ch`. Where `ch` is
[`ConvergenceHistory`](@ref) object. Otherwise it will only return `x`.

The `plot` attribute can only be used when using the `Master` version.

**Arguments**

* `K::KrylovSubspace`: krylov subspace.

* `A`: linear operator.

* `Master::Type{Master}`: dispatch type.

*Keywords*

$karg

* `x = random unit vector`: initial eigenvector guess.

* `tol::Real = eps()*size(A,2)^3`: stopping tolerance.

* `maxiter::Integer = size(A,2)`: maximum number of iterations.

* `verbose::Bool = false`: verbose flag.

* `plot::Bool = false`: plot data. (Only with `Master` version)

**Output**

* `::Real`: eigen value

* `::Vector`: eigen vector

*ConvergenceHistory keys*

* `:tol` => `::Real`: stopping tolerance.

* `:resnom` => `::Vector`: residual norm at each iteration.

""" -> func
end
