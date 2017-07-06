export cg, cg!

cg(A, b; kwargs...) = cg!(zerox(A, b), A, b; kwargs...)

function cg!(x, A, b;
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Integer = min(20, size(A, 1)),
    plot = false,
    log::Bool = false,
    Pl = Identity(),
    kwargs...
)
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && reserve!(history, :resnorm, maxiter + 1)
    cg_method!(history, x, A, b, Pl; tol = tol, log = log, maxiter = maxiter, kwargs...)
    log && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

function cg_method!(history::ConvergenceHistory, x, A, b, Pl;
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Integer = min(20, size(A, 1)),
    verbose::Bool = false,
    log = false
)
    # Preconditioned CG
    T = eltype(b)
    n = size(A, 1)

    # Initial residual vector
    r = copy(b)
    c = A * x
    @blas! r -= one(T) * c
    u = zeros(T, n)
    ρ = one(T)

    if log
        history.mvps += 1
    end
    
    iter = 1

    # Here you could save one inner product if norm(r) is used rather than norm(b)
    reltol = norm(b) * tol
    last_residual = zero(T)

    while true

        last_residual = norm(r)

        verbose && @printf("%3d\t%1.2e\n", iter, last_residual)

        if last_residual ≤ reltol || iter > maxiter
            break
        end

        # Log progress        
        if log
            nextiter!(history, mvps = 1)
            push!(history, :resnorm, last_residual)
        end

        # Preconditioner: c = Pl \ r
        solve!(c, Pl, r)

        ρ_prev = ρ
        ρ = dot(c, r)
        β = ρ / ρ_prev

        # u := c + βu (almost an axpy)
        @blas! u *= β
        @blas! u += one(T) * c

        # c = A * u
        A_mul_B!(c, A, u)
        α = ρ / dot(u, c)
    
        # Improve solution and residual
        @blas! x += α * u
        @blas! r -= α * c

        iter += 1
    end

    verbose && @printf("\n")
    log && setconv(history, last_residual < reltol)

    x
end

function cg_method!(history::ConvergenceHistory, x, A, b, Pl::Identity;
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Integer = min(20, size(A, 1)),
    verbose::Bool = false,
    log = false
)
    # Unpreconditioned CG
    T = eltype(b)
    n = size(A, 1)

    # Initial residual vector
    r = copy(b)
    c = A * x
    @blas! r -= one(T) * c
    u = zeros(T, n)
    ρ = one(T)

    if log
        history.mvps += 1
    end

    iter = 1

    reltol = norm(b) * tol
    last_residual = zero(T)

    while true

        ρ_prev = ρ
        ρ = dot(r, r)
        β = ρ / ρ_prev

        last_residual = sqrt(ρ)

        # Log progress
        if log
            nextiter!(history, mvps = 1)
            push!(history, :resnorm, last_residual)
        end

        verbose && @printf("%3d\t%1.2e\n", iter, last_residual)

        # Stopping condition
        if last_residual ≤ reltol || iter > maxiter
            break
        end

        # u := r + βu (almost an axpy)
        @blas! u *= β
        @blas! u += one(T) * r

        # c = A * u
        A_mul_B!(c, A, u)
        α = ρ / dot(u, c)
    
        # Improve solution and residual
        @blas! x += α * u
        @blas! r -= α * c

        iter += 1
    end

    verbose && @printf("\n")
    log && setconv(history, last_residual < reltol)

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

`Pl = Identity()`: left preconditioner of the method.

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