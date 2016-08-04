export cg, cg!

####################
# API method calls #
####################

cg(A, b; kwargs...) =  cg!(zerox(A,b), A, b; kwargs...)
cg(::Type{Master}, A, b; kwargs...) =  cg!(Master, zerox(A,b), A, b; kwargs...)


function cg!(x, A, b; kwargs...)
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(x,K,b; kwargs...)
end
function cg!(x, K::KrylovSubspace, b; kwargs...)
    cg_method!(x,K,b; kwargs...)
    x
end
function cg!(::Type{Master}, x, A, b; kwargs...)
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(Master, x,K,b; kwargs...)
end
function cg!(::Type{Master}, x, K::KrylovSubspace, b;
    tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2),
    verbose::Bool=false, plot=false, Pl=1
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm, maxiter)
    cg_method!(x,K,b; Pl=Pl,tol=tol,maxiter=maxiter,log=log,verbose=verbose)
    shrink!(log)
    plot && showplot(log)
    x, log
end

#########################
# Method Implementation #
#########################

function cg_method!(x,K,b;
    Pl=1,tol::Real=size(K.A,2)*eps(),maxiter::Integer=size(K.A,2),
    log::MethodLog=DummyHistory(),verbose::Bool=false
    )
    verbose && @printf("=== cg ===\n%4s\t%7s\n","iter","relres")
    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = Pl\r
    γ = dot(r, z)
    for iter=1:maxiter
        append!(K, p)
        q = nextvec(K)
        α = γ/dot(p, q)
        # α>=0 || throw(PosSemidefException("α=$α"))
        update!(x, α, p)
        r -= α*q
        resnorm = norm(r)
        nextiter!(log)
        push!(log,:resnorm,resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
        z = Pl\r
        oldγ = γ
        γ = dot(r, z)
        β = γ/oldγ
        p = z + β*p
    end
    setmvps(log, K.mvps)
    setconv(log, 0<=norm(r)<tol)
    verbose && @printf("\n")
end

#################
# Documentation #
#################

#Initialize parameters
doc_call = """    cg(A, b)
    cg(Master, A, b)
"""
doc!_call = """    cg!(x, A, b)
    cg!(x, K, b)
    cg!(Master, x, A, b)
    cg!(Master, x, K, b)
"""

doc_msg = "Solve A*x=b with the conjugate gradients method."
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """* `x`: initial guess, overwrite final estimation.

* `K::KrylovSubspace`: krylov subspace."""

doc_version = (cg, doc_call, doc_msg, doc_arg)
doc!_version = (cg!, doc!_call, doc!_msg, doc!_arg)

#Build docs
for (func, call, msg, arg) in [doc_version, doc!_version]
@doc """
$call

$msg

If [`Master`](@ref) is given, method will output a tuple `x, ch`. Where `ch` is
[`ConvergenceHistory`](@ref) object. Otherwise it will only return `x`.

The `plot` attribute can only be used when using the `Master` version.

**Arguments**

$arg

* `A`: linear operator.

* `b`: right hand side.

* `Master::Type{Master}`: dispatch type.

*Keywords*

* `Pl = 1`: left preconditioner of the method.

* `tol::Real = size(A,2)*eps()`: stopping tolerance.

* `maxiter::Integer = size(A,2)`: maximum number of iterations.

* `verbose::Bool = false`: print method information.

* `plot::Bool = false`: plot data. (Only with `Master` version)

**Output**

*Normal version:*

* `x`: approximated solution.

*`Master` version:*

* `x`: approximated solution.

* `ch`: convergence history.

*ConvergenceHistory keys*

* `:tol` => `::Real`: stopping tolerance.

* `:resnom` => `::Vector`: residual norm at each iteration.

""" -> func
end
