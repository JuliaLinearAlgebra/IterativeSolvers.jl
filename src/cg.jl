export cg, cg!

####################
# API method calls #
####################

cg(A, b; kwargs...) = cg!(zerox(A,b), A, b; kwargs...)

function cg!(x, A, b;
    tol::Real=size(A,2)*eps(), maxiter::Integer=size(A,2),
    plot=false, log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    cg_method!(history, x, K, b; tol=tol, maxiter=maxiter, kwargs...)
    plot && (shrink!(history); showplot(history))
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

function cg_method!(log::ConvergenceHistory, x, K, b;
    Pl=1,tol::Real=size(K.A,2)*eps(),maxiter::Integer=size(K.A,2), verbose::Bool=false
    )
    verbose && @printf("=== cg ===\n%4s\t%7s\n","iter","resnorm")
    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = Pl\r
    γ = dot(r, z)
    for iter=1:maxiter
        nextiter!(log, mvps=1)
        append!(K, p)
        q = nextvec(K)
        α = γ/dot(p, q)
        # α>=0 || throw(PosSemidefException("α=$α"))
        update!(x, α, p)
        r -= α*q
        resnorm = norm(r)
        push!(log,:resnorm,resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
        z = Pl\r
        oldγ = γ
        γ = dot(r, z)
        β = γ/oldγ
        p = z + β*p
    end
    verbose && @printf("\n")
    setconv(log, 0<=norm(r)<tol)
    x
end

#################
# Documentation #
#################

let
#Initialize parameters
doc_call = """    cg(A, b)
"""
doc!_call = """    cg!(x, A, b)
"""

doc_msg = "Solve A*x=b with the conjugate gradients method."
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """* `x`: initial guess, overwrite final estimation."""

doc_version = (doc_call, doc_msg, doc_arg)
doc!_version = (doc!_call, doc!_msg, doc!_arg)

i = 0
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

* `Pl = 1`: left preconditioner of the method.
* `tol::Real = size(A,2)*eps()`: stopping tolerance.
* `maxiter::Integer = size(A,2)`: maximum number of iterations.
* `verbose::Bool = false`: print method information.
* `log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.
* `plot::Bool = false`: plot data. (Only when `log` is set)

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

@doc docstring[1] -> cg
@doc docstring[2] -> cg!
end
