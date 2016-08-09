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
function powm{T}(K::KrylovSubspace{T};
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false, log::Bool=false
    )
    if log
        history = ConvergenceHistory()
        history[:tol] = tol
        reserve!(history,:resnorm, maxiter)
    else
        history = DummyHistory()
    end
    eig, v = powm_method(K; verbose=verbose,log=history,tol=tol,maxiter=maxiter)
    if log
        shrink!(history)
        plot && showplot(history)
        eig, v, history
    else
        eig, v
    end
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
function invpowm{T}(K::KrylovSubspace{T};
    shift::Number=0, tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false, log::Bool=false
    )
    if log
        history = ConvergenceHistory()
        history[:tol] = tol
        reserve!(history,:resnorm, maxiter)
    else
        history = DummyHistory()
    end
    eig, v = invpowm_method(K, shift;
        verbose=verbose, log=history, tol=tol, maxiter=maxiter
        )
    if log
        shrink!(history)
        plot && showplot(history)
        eig, v, history
    else
        eig, v
    end
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
function rqi{T}(K::KrylovSubspace{T};
    shift::Number=0, tol::Real=eps(T)*K.n^3, maxiter::Int=K.n,
    plot::Bool=false, verbose::Bool=false, log::Bool=false
    )
    if log
        history = ConvergenceHistory()
        history[:tol] = tol
        reserve!(history,:resnorm, maxiter)
    else
        history = DummyHistory()
    end
    eig, v = rqi_method(K, shift; verbose=verbose, log=history, tol=tol, maxiter=maxiter)
    if log
        shrink!(history)
        plot && showplot(history)
        eig, v, history
    else
        eig, v
    end
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
"""
doc2_call = """    invpowm(A)
    invpowm(K)
"""
doc3_call = """    rqi(A)
    rqi(K)
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

If `log` is set to `true` is given, method will output a tuple `eig, v, ch`. Where
`ch` is a [`ConvergenceHistory`](@ref) object. Otherwise it will only return `eig, v`.

The `plot` attribute can only be used when `log` is set version.

**Arguments**

* `K::KrylovSubspace`: krylov subspace.

* `A`: linear operator.

*Keywords*

$karg

* `x = random unit vector`: initial eigenvector guess.

* `tol::Real = eps()*size(A,2)^3`: stopping tolerance.

* `maxiter::Integer = size(A,2)`: maximum number of iterations.

* `verbose::Bool = false`: verbose flag.

* `log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

* `plot::Bool = false`: plot data. (Only when `log` is set)

**Output**

*`log` is `false`:*

* `eig::Real`: eigen value

* `v::Vector`: eigen vector

*`log` is `true`:*

* `eig::Real`: eigen value

* `v::Vector`: eigen vector

* `ch`: convergence history.

*ConvergenceHistory keys*

* `:tol` => `::Real`: stopping tolerance.

* `:resnom` => `::Vector`: residual norm at each iteration.

""" -> func
end
