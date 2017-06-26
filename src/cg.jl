export cg, cg!

cg(A, b; kwargs...) = cg!(zerox(A, b), A, b; kwargs...)

function cg!(x, A, b;
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Integer = min(20, size(A, 1)),
    plot = false,
    log::Bool = false,
    kwargs...
)
    (plot & !log) && error("Can't plot when log keyword is false")
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    reserve!(history, :resnorm, maxiter)
    cg_method!(history, x, K, b; tol = tol, maxiter = maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

function cg_method!(log::ConvergenceHistory, x, K, b;
    Pl = 1,
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Integer = min(20, size(K.A, 1)),
    verbose::Bool = false
)
    T = eltype(b)
    n = size(K.A, 1)

    # Initial residual vector
    r = copy(b)
    @blas! r -= one(T) * nextvec(K)
    c = zeros(T, n)
    u = zeros(T, n)
    ρ = one(T)
    
    iter = 0

    last_residual = norm(r)

    # Here you could save one inner product if norm(r) is used rather than norm(b)
    reltol = norm(b) * tol

    while last_residual > reltol && iter < maxiter
        nextiter!(log, mvps = 1)

        # Preconditioner: c = Pl \ r
        solve!(c, Pl, r)

        ρ_prev = ρ
        ρ = dot(c, r)
        β = -ρ / ρ_prev

        # u := r - βu (almost an axpy)
        @blas! u *= -β
        @blas! u += one(T) * c

        # c = A * u
        append!(K, u)
        nextvec!(c, K)
        α = ρ / dot(u, c)
    
        # Improve solution and residual
        @blas! x += α * u
        @blas! r -= α * c

        iter += 1
        last_residual = norm(r)

        # Log progress
        push!(log, :resnorm, last_residual)
        verbose && @printf("%3d\t%1.2e\n", iter, last_residual)
    end

    verbose && @printf("\n")
    setconv(log, last_residual < reltol)

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
doc!_arg = """`x`: initial guess, overwrite final estimation."""

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

# Arguments

$arg

`A`: linear operator.

`b`: right hand side.

## Keywords

`Pl = 1`: left preconditioner of the method.

`tol::Real = size(A,2)*eps()`: stopping tolerance.

`maxiter::Integer = size(A,2)`: maximum number of iterations.

`verbose::Bool = false`: print method information.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

`plot::Bool = false`: plot data. (Only when `log` is set)

# Output

**if `log` is `false`**

`x`: approximated solution.

**if `log` is `true`**

`x`: approximated solution.

`ch`: convergence history.

**ConvergenceHistory keys**

`:tol` => `::Real`: stopping tolerance.

`:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> cg
@doc docstring[2] -> cg!
end