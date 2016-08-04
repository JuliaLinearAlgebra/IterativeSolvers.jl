import Base.LinAlg.BlasFloat

export eiglancz

####################
# API method calls #
####################

eiglancz(A; kwargs...) = eiglancz_method(A; kwargs...)
function eiglancz(::Type{Master}, A;
    maxiter::Integer=size(A,1), plot::Bool=false,
    tol::Real = size(A,1)^3*eps(real(eltype(A))), kwargs...
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log, :resnorm,maxiter)
    eigs = eiglancz_method(A; maxiter=maxiter, log=log, kwargs...)
    shrink!(log)
    plot && showplot(log)
    eigs, log
end

#########################
# Method Implementation #
#########################

function lanczos!{T}(K::KrylovSubspace{T})
    m = K.n
    αs = Array(T, m)
    βs = Array(T, m-1)
    for j=1:m-1
        w = nextvec(K)
        if j>1 w -= βs[j-1]*K.v[1] end
        w, y = orthogonalize(w, K, 1)
        αs[j] = y[1]
        βs[j] = convert(T, norm(w))
        append!(K, w/βs[j])
    end
    αs[m]= dot(nextvec(K), lastvec(K))
    SymTridiagonal(αs, βs)
end

function eiglancz_method(A;
    neigs::Int=size(A,1), tol::Real = size(A,1)^3*eps(real(eltype(A))),
    maxiter::Integer=size(A,1), verbose::Bool=false, log::MethodLog=DummyHistory()
    )
    verbose && @printf("=== eiglancz ===\n%4s\t%7s\n","iter","relres")
    K = KrylovSubspace(A, size(A, 1), 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    e1 = eigvals(lanczos!(K), 1:neigs)
    for iter=1:maxiter
        e0, e1 = e1, eigvals(lanczos!(K), 1:neigs)
        resnorm = norm(e1-e0)
        nextiter!(log)
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
    end
    setmvps(log, K.mvps)
    verbose && @printf("\n")
    e1
end

#################
# Documentation #
#################

#Initialize parameters
doc_call = """    eiglancz(A)
"""

doc_msg = "Find the most useful eigenvalues using the lanczos method."

doc_arg = ""

doc_version = (eiglancz, doc_call, doc_msg, doc_arg)

#Build docs
for (func, call, msg, arg) in [doc_version]
@doc """
$call

$msg

If [`Master`](@ref) is given, method will output a tuple `x, ch`. Where `ch` is
[`ConvergenceHistory`](@ref) object. Otherwise it will only return `x`.

The `plot` attribute can only be used when using the `Master` version.

**Arguments**

$arg

* `A`: linear operator.

* `Master::Type{Master}`: dispatch type.

*Keywords*

* `neigs::Int = size(A,1)`: number of eigen values.

* `tol::Real = size(A,1)^3*eps()`: stopping tolerance.

* `maxiter::Integer=size(A,1)`: maximum number of iterations.

* `verbose::Bool = false`: verbose flag.

* `plot::Bool = false`: plot data. (Only with `Master` version)

**Output**

* `::Vector`: eigen values.

*ConvergenceHistory keys*

* `:tol` => `::Real`: stopping tolerance.

* `:resnom` => `::Vector`: residual norm at each iteration.

""" -> func
end
