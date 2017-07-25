import Base: next, start, done

export chebyshev, chebyshev!

type ChebyshevIterable{precT, matT, vecT, realT <: Real}
    Pl::precT
    A::matT
    b::vecT

    x::vecT
    r::vecT
    u::vecT
    c::vecT

    α::realT

    λ_avg::realT
    λ_diff::realT

    resnorm::realT
    reltol::realT
    maxiter::Int
    mv_products::Int
end

converged(c::ChebyshevIterable) = c.resnorm ≤ c.reltol
start(::ChebyshevIterable) = 0
done(c::ChebyshevIterable, iteration::Int) = iteration ≥ c.maxiter || converged(c)

function next(cheb::ChebyshevIterable, iteration::Int)
    T = eltype(cheb.u)

    solve!(cheb.c, cheb.Pl, cheb.r)

    if iteration == 1
        cheb.α = T(2) / cheb.λ_avg
        copy!(cheb.u, cheb.c)
    else
        β = (cheb.λ_diff * cheb.α / 2) ^ 2
        cheb.α = inv(cheb.λ_avg - β)

        # Almost an axpy u = c + βu
        scale!(cheb.u, β)
        @blas! cheb.u += one(T) * cheb.c
    end

    A_mul_B!(cheb.c, cheb.A, cheb.u)
    cheb.mv_products += 1

    @blas! cheb.x += cheb.α * cheb.u
    @blas! cheb.r -= cheb.α * cheb.c

    cheb.resnorm = norm(cheb.r)

    cheb.resnorm, iteration + 1
end

chebyshev_iterable(A, b, λmin::Real, λmax::Real; kwargs...) =
    chebyshev_iterable!(zerox(A, b), A, b, λmin, λmax; kwargs...)

function chebyshev_iterable!(x, A, b, λmin::Real, λmax::Real;
    tol = sqrt(eps(real(eltype(b)))), maxiter = size(A, 1), Pl = Identity(), initially_zero = false)

    λ_avg = (λmax + λmin) / 2
    λ_diff = (λmax - λmin) / 2

    T = eltype(b)
    r = copy(b)
    u = zeros(x)
    c = similar(x)

    # One MV product less
    if initially_zero
        resnorm = norm(r)
        reltol = tol * resnorm
        mv_products = 0
    else
        A_mul_B!(c, A, x)
        @blas! r -= one(T) * c
        resnorm = norm(r)
        reltol = tol * norm(b)
        mv_products = 1
    end

    ChebyshevIterable(Pl, A, b,
        x, r, u, c,
        zero(real(T)),
        λ_avg, λ_diff,
        resnorm, reltol, maxiter, mv_products
    )
end

####################
# API method calls #
####################

chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) =
    chebyshev!(zerox(A, b), A, b, λmin, λmax; initially_zero = true, kwargs...)

function chebyshev!(x, A, b, λmin::Real, λmax::Real;
    Pl = Identity(),
    tol::Real=sqrt(eps(real(eltype(b)))),
    maxiter::Int=size(A, 1),
    log::Bool=false,
    verbose::Bool=false,
    initially_zero::Bool=false
)
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history, :resnorm, maxiter)

    verbose && @printf("=== chebyshev ===\n%4s\t%7s\n","iter","resnorm")

    iterable = chebyshev_iterable!(x, A, b, λmin, λmax; tol=tol, maxiter=maxiter, Pl=Pl, initially_zero=initially_zero)
    history.mvps = iterable.mv_products
    for (iteration, resnorm) = enumerate(iterable)
        nextiter!(history)
        history.mvps = iterable.mv_products
        push!(history, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n", iteration, resnorm)
    end
    verbose && println()
    setconv(history, converged(iterable))
    log && shrink!(history)

    log ? (x, history) : x
end

#################
# Documentation #
#################

let
#Initialize parameters
doc_call = """    chebyshev!(A, b, λmin, λmax)
"""
doc!_call = """    chebyshev!(x, A, b, λmin, λmax)
"""

doc_msg = "Solve A*x=b using the chebyshev method."
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """`x`: initial guess, overwrite final estimation."""

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

# Arguments

$arg

`A`: linear operator.

`b`: right hand side.

## Keywords

`Pl = 1`: left preconditioner of the method.

`tol::Real = sqrt(eps())`: stopping tolerance.

`maxiter::Integer = size(A,2)^3`: maximum number of iterations.

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

@doc docstring[1] -> chebyshev
@doc docstring[2] -> chebyshev!
end
