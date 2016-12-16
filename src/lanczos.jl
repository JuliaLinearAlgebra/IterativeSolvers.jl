import Base.LinAlg.BlasFloat

export eiglancz

####################
# API method calls #
####################

function eiglancz(A;
    maxiter::Integer=size(A,1), plot::Bool=false,
    tol::Real = size(A,1)^3*eps(real(eltype(A))), log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    e1 = eiglancz_method(history, A; maxiter=maxiter, tol=tol, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (e1, history) : e1
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

function eiglancz_method(log::ConvergenceHistory, A;
    neigs::Int=size(A,1), tol::Real = size(A,1)^3*eps(real(eltype(A))),
    maxiter::Integer=size(A,1), verbose::Bool=false
    )
    verbose && @printf("=== eiglancz ===\n%4s\t%7s\n","iter","resnorm")
    K = KrylovSubspace(A, size(A, 1), 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    e1 = eigvals(lanczos!(K), 1:neigs)
    for iter=1:maxiter
        nextiter!(log,mvps=1)
        e0, e1 = e1, eigvals(lanczos!(K), 1:neigs)
        resnorm = norm(e1-e0)
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
    end
    verbose && @printf("\n")
    e1
end

#################
# Documentation #
#################

let
#Initialize parameters
doc_call = """    eiglancz(A)
"""

doc_msg = "Find the most useful eigenvalues using the lanczos method."

doc_arg = ""

doc_version = (eiglancz, doc_call, doc_msg, doc_arg)

i=0
docstring = Vector(1)

#Build docs
for (func, call, msg, arg) in [doc_version]
i+=1
docstring[i] = """
$call

$msg

If `log` is set to `true` is given, method will output a tuple `eigs, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `eigs`.

The `plot` attribute can only be used when `log` is set version.

# Arguments

$arg

`A`: linear operator.

## Keywords

`neigs::Int = size(A,1)`: number of eigen values.

`tol::Real = size(A,1)^3*eps()`: stopping tolerance.

`maxiter::Integer=size(A,1)`: maximum number of iterations.

`verbose::Bool = false`: verbose flag.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

`plot::Bool = false`: plot data. (Only when `log` is set)

# Output

**if `log` is `false`**

`eigs::Vector`: eigen values.

**if `log` is `true`**

`eigs::Vector`: eigen values.

`ch`: convergence history.

**ConvergenceHistory keys**

`:tol` => `::Real`: stopping tolerance.

`:resnom` => `::Vector`: residual norm at each iteration.
"""
end

@doc docstring[1] -> eiglancz
end
