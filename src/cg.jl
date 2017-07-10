import Base: start, next, done

mutable struct CGIterable{matT <: AbstractMatrix, vecT <: AbstractVector, numT <: Real}
    A::matT
    x::vecT
    b::vecT
    r::vecT
    c::vecT
    u::vecT
    reltol::numT
    residual::numT
    prev_residual::numT
    maxiter::Int
    mv_products::Int
end

mutable struct PCGIterable{precT, matT <: AbstractMatrix, vecT <: AbstractVector, numT <: Real, paramT <: Number}
    Pl::precT
    A::matT
    x::vecT
    b::vecT
    r::vecT
    c::vecT
    u::vecT
    reltol::numT
    residual::numT
    ρ::paramT
    maxiter::Int
    mv_products::Int
end

@inline converged(it::Union{CGIterable, PCGIterable}) = it.residual ≤ it.reltol

@inline start(it::Union{CGIterable, PCGIterable}) = 0

@inline done(it::Union{CGIterable, PCGIterable}, iteration::Int) = iteration ≥ it.maxiter || converged(it)


###############
# Ordinary CG #
###############

function next(it::CGIterable, iteration::Int)
    # u := r + βu (almost an axpy)
    β = it.residual^2 / it.prev_residual^2
    @blas! it.u *= β
    @blas! it.u += one(eltype(it.b)) * it.r

    # c = A * u
    A_mul_B!(it.c, it.A, it.u)
    α = it.residual^2 / dot(it.u, it.c)

    # Improve solution and residual
    @blas! it.x += α * it.u
    @blas! it.r -= α * it.c

    it.prev_residual = it.residual
    it.residual = norm(it.r)

    # Return the residual at item and iteration number as state
    it.residual, iteration + 1
end

#####################
# Preconditioned CG #
#####################

function next(it::PCGIterable, iteration::Int)
    A_ldiv_B!(it.c, it.Pl, it.r)

    ρ_prev = it.ρ
    it.ρ = dot(it.c, it.r)
    
    # u := c + βu (almost an axpy)
    β = it.ρ / ρ_prev
    @blas! it.u *= β
    @blas! it.u += one(eltype(it.b)) * it.c

    # c = A * u
    A_mul_B!(it.c, it.A, it.u)
    α = it.ρ / dot(it.u, it.c)

    # Improve solution and residual
    @blas! it.x += α * it.u
    @blas! it.r -= α * it.c

    it.residual = norm(it.r)

    # Return the residual at item and iteration number as state
    it.residual, iteration + 1
end

# Utility functions

@inline cg_iterator(A, b, Pl = Identity(); kwargs...) = cg_iterator!(zerox(A, b), A, b, Pl; initially_zero = true, kwargs...)

function cg_iterator!(x, A, b, Pl = Identity();
    tol = sqrt(eps(real(eltype(b)))),
    maxiter = min(20, length(b)),
    initially_zero::Bool = false
)
    u = zeros(x)
    r = copy(b)

    # Compute r with an MV-product or not.
    if initially_zero
        mv_products = 0
        c = similar(x)
    else
        mv_products = 1
        c = A * x
        @blas! r -= one(eltype(x)) * c
    end

    # Stopping criterion
    residual = norm(r)
    ρ = one(residual)
    reltol = residual * tol

    # Return the iterable
    if isa(Pl, Identity)
        return CGIterable(A, x, b,
            r, c, u,
            reltol, residual, ρ,
            maxiter, mv_products
        )
    else
        return PCGIterable(Pl, A, x, b,
            r, c, u,
            reltol, residual, ρ,
            maxiter, mv_products
        )
    end
end

cg(A, b, Pl = Identity(); kwargs...) = cg!(zerox(A, b), A, b, Pl; initially_zero = true, kwargs...)

function cg!(x, A, b;
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Integer = min(20, size(A, 1)),
    plot = false,
    log::Bool = false,
    verbose::Bool = false,
    Pl = Identity(),
    kwargs...
)
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && reserve!(history, :resnorm, maxiter + 1)

    # Actually perform CG
    iterable = cg_iterator!(x, A, b, Pl; tol = tol, maxiter = maxiter, kwargs...)
    if log
        history.mvps = iterable.mv_products
    end
    for (iteration, item) = enumerate(iterable)
        if log
            nextiter!(history, mvps = 1)
            push!(history, :resnorm, iterable.residual)
        end
        verbose && @printf("%3d\t%1.2e\n", iteration, iterable.residual)
    end

    verbose && println()
    log && setconv(history, converged(iterable))
    log && shrink!(history)
    plot && showplot(history)

    log ? (iterable.x, history) : iterable.x
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