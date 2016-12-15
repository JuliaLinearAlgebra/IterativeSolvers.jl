#Simple methods
export powm, invpowm, rqi

####################
# API method calls #
####################

function powm(A;
    x=nothing, tol::Real=eps(eltype(A))*size(A,2)^3, maxiter::Int=size(A,2),
    plot::Bool=false, log::Bool=false, kwargs...
    )
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))

    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    eig, v = powm_method!(history, K; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (eig, v, history) : (eig, v)
end

#########################
# Method Implementation #
#########################

function powm_method!{T}(log::ConvergenceHistory, K::KrylovSubspace{T};
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n, verbose::Bool=false
    )
    verbose && @printf("=== powm ===\n%4s\t%7s\n","iter","resnorm")
    θ = zero(T)
    v = Array(T, K.n)
    for iter=1:maxiter
        nextiter!(log,mvps=1)
        v = lastvec(K)
        y = nextvec(K)
        θ = dot(v, y)
        resnorm = real(norm(y - θ*v))
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm <= tol*abs(θ) && (setconv(log, resnorm >= 0); break)
        appendunit!(K, y)
    end
    verbose && @printf("\n")
    θ, v
end

####################
# API method calls #
####################

function invpowm(A;
    x=nothing, shift::Number=0, tol::Real=eps(eltype(A))*size(A,2)^3,
    maxiter::Int=size(A,2), plot::Bool=false, log::Bool=false, kwargs...
    )
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))

    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    eig, v = invpowm_method!(history, K, shift; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (eig, v, history) : (eig, v)
end

#########################
# Method Implementation #
#########################

function invpowm_method!{T}(log::ConvergenceHistory, K::KrylovSubspace{T}, σ::Number=zero(T);
    tol::Real=eps(T)*K.n^3, maxiter::Int=K.n, verbose::Bool=false
    )
    verbose && @printf("=== invpowm ===\n%4s\t%7s\n","iter","resnorm")
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
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm <= tol*abs(θ) && (setconv(log, resnorm >= 0); break)
        appendunit!(K, y)
    end
    verbose && @printf("\n")
    σ+1/θ, y/θ
end

####################
# API method calls #
####################

function rqi(A;
    x=nothing, shift::Number=0, tol::Real=eps(eltype(A))*size(A,2)^3,
    maxiter::Int=size(A,2), plot::Bool=false, log::Bool=false, kwargs...
    )
    K = KrylovSubspace(A, 1)
    x==nothing ? initrand!(K) : init!(K, x/norm(x))

    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    eig, v = rqi_method!(history, K, shift; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (eig, v, history) : (eig, v)
end

#########################
# Method Implementation #
#########################

#Rayleigh quotient iteration
#XXX Doesn't work well
function rqi_method!{T}(log::ConvergenceHistory, K::KrylovSubspace{T}, σ::Number;
    tol::Real=eps(T), maxiter::Int=K.n, verbose::Bool=false
    )
    verbose && @printf("=== rqi ===\n%4s\t%7s\n","iter","resnorm")
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
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        θ >= 1/tol && (setconv(log, resnorm >= 0); break)
    end
    verbose && @printf("\n")
    ρ, v
end

#################
# Documentation #
#################

let
#Initialize parameters
doc1_call = """    powm(A)
"""
doc2_call = """    invpowm(A)
"""
doc3_call = """    rqi(A)
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
doc2_karg = "`shift::Number=0`: shift to be applied to matrix A."
doc3_karg = "`shift::Number=0`: shift to be applied to matrix A."

doc1_version = (powm, doc1_call, doc1_msg, doc1_karg)
doc2_version = (invpowm, doc2_call, doc2_msg, doc2_karg)
doc3_version = (rqi, doc3_call, doc3_msg, doc3_karg)

i=0
docstring = Vector(3)

#Build docs
for (func, call, msg, karg) in [doc1_version, doc2_version, doc3_version]
i+=1
docstring[i] = """
$call

$msg

If `log` is set to `true` is given, method will output a tuple `eig, v, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `eig, v`.

The `plot` attribute can only be used when `log` is set version.

# Arguments

`K::KrylovSubspace`: krylov subspace.

`A`: linear operator.

## Keywords

$karg

`x = random unit vector`: initial eigenvector guess.

`tol::Real = eps()*size(A,2)^3`: stopping tolerance.

`maxiter::Integer = size(A,2)`: maximum number of iterations.

`verbose::Bool = false`: verbose flag.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

`plot::Bool = false`: plot data. (Only when `log` is set)

# Output

**if `log` is `false`**

`eig::Real`: eigen value

`v::Vector`: eigen vector

**if `log` is `true`**

`eig::Real`: eigen value

`v::Vector`: eigen vector

`ch`: convergence history.

**ConvergenceHistory keys**

`:tol` => `::Real`: stopping tolerance.

`:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> powm
@doc docstring[2] -> invpowm
@doc docstring[3] -> rqi
end
