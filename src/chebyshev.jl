export chebyshev, chebyshev!

####################
# API method calls #
####################

chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) =
    chebyshev!(zerox(A, b), A, b, λmin, λmax; kwargs...)

function chebyshev!(x, A, b, λmin::Real, λmax::Real; n::Int=size(A,2), kwargs...)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	chebyshev!(x, K, b, λmin, λmax; kwargs...)
end
function chebyshev!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
    n::Int=size(K,2), tol::Real = sqrt(eps(typeof(real(b[1])))),
    maxiter::Int = n^3, plot::Bool=false, log::Bool=false, kwargs...
    )
    if log
        history = ConvergenceHistory()
        history[:tol] = tol
        reserve!(history,:resnorm,maxiter)
    else
        history = DummyHistory()
    end
    chebyshev_method!(x,K,b,λmin,λmax; tol=tol,maxiter=maxiter,log=history,kwargs...)
    if log
        shrink!(history)
        plot && showplot(history)
        x, history
    else
        x
    end
end

#########################
# Method Implementation #
#########################

function chebyshev_method!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
    Pr=1, tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = K.n^3,
    log::MethodLog=DummyHistory(), verbose::Bool=false
    )
    verbose && @printf("=== chebyshev ===\n%4s\t%7s\n","iter","relres")
    local α, p
    K.order = 1
    tol = tol*norm(b)
    r = b - nextvec(K)
    d::eltype(b) = (λmax + λmin)/2
    c::eltype(b) = (λmax - λmin)/2
    for iter = 1:maxiter
    	z = Pr\r
    	if iter == 1
    		p = z
    		α = 2/d
    	else
    		β = (c*α/2)^2
    		α = 1/(d - β)
    		p = z + β*p
    	end
    	append!(K, p)
    	update!(x, α, p)
    	r -= α*nextvec(K)
        resnorm = norm(r)
        nextiter!(log)
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
	end
    setmvps(log, K.mvps)
    setconv(log, 0<=norm(r)<tol)
    verbose && @printf("\n")
end

#################
# Documentation #
#################

#Initialize parameters
doc_call = """    chebyshev!(A, b, λmin, λmax)
"""
doc!_call = """    chebyshev!(x, A, b, λmin, λmax)
    chebyshev!(x, K, b, λmin, λmax)
"""

doc_msg = "Solve A*x=b using the chebyshev method."
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """* `x`: initial guess, overwrite final estimation.

* `K::KrylovSubspace`: krylov subspace."""

doc_version = (doc_call, doc_msg, doc_arg)
doc!_version = (doc!_call, doc!_msg, doc!_arg)

i=0
docstring = Vector(2)

#Build docs
for (call, msg, arg) in [doc_version, doc!_version] #Start
i+=1
docstring[i] = """
$call

$msg

If `log` is set to `true` is given, method will output a tuple `x, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `x`.

The `plot` attribute can only be used when `log` is set version.

**Arguments**

$arg

* `A`: linear operator.

* `b`: right hand side.

*Keywords*

* `Pr = 1`: right preconditioner of the method.

* `tol::Real = sqrt(eps())`: stopping tolerance.

* `maxiter::Integer = size(A,2)^3`: maximum number of iterations.

* `verbose::Bool = false`: print method information.

* `log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

* `plot::Bool = false`: plot data. (Only with `Master` version)

**Output**

*`log` is `false`:*

* `x`: approximated solution.

*`log` is `true`:*

* `x`: approximated solution.

* `ch`: convergence history.

*ConvergenceHistory keys*

* `:tol` => `::Real`: stopping tolerance.

* `:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> chebyshev
@doc docstring[2] -> chebyshev!
