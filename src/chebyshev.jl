import Base: next, start, done

export chebyshev, chebyshev!

mutable struct ChebyshevIterable{precT, matT, solT, vecT, realT <: Real}
    Pl::precT
    A::matT

    x::solT
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
    T = eltype(cheb.x)

    A_ldiv_B!(cheb.c, cheb.Pl, cheb.r)

    if iteration == 1
        cheb.α = T(2) / cheb.λ_avg
        copyto!(cheb.u, cheb.c)
    else
        β = (cheb.λ_diff * cheb.α / 2) ^ 2
        cheb.α = inv(cheb.λ_avg - β)
        cheb.u .= cheb.c .+ β .* cheb.c
    end

    A_mul_B!(cheb.c, cheb.A, cheb.u)
    cheb.mv_products += 1

    cheb.x .+= cheb.α .* cheb.u
    cheb.r .-= cheb.α .* cheb.c

    cheb.resnorm = norm(cheb.r)

    cheb.resnorm, iteration + 1
end

function chebyshev_iterable!(x, A, b, λmin::Real, λmax::Real;
    tol = sqrt(eps(real(eltype(b)))), maxiter = size(A, 2), Pl = Identity(), initially_zero = false)

    λ_avg = (λmax + λmin) / 2
    λ_diff = (λmax - λmin) / 2

    T = eltype(x)
    r = similar(x)
    copyto!(r, b)
    u = zeros(x)
    c = similar(x)

    # One MV product less
    if initially_zero
        resnorm = norm(r)
        reltol = tol * resnorm
        mv_products = 0
    else
        A_mul_B!(c, A, x)
        r .-= c
        resnorm = norm(r)
        reltol = tol * norm(b)
        mv_products = 1
    end

    ChebyshevIterable(Pl, A, x, r, u, c,
        zero(real(T)),
        λ_avg, λ_diff,
        resnorm, reltol, maxiter, mv_products
    )
end


"""
    chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]

Same as [`chebyshev!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) =
    chebyshev!(zerox(A, b), A, b, λmin, λmax; initially_zero = true, kwargs...)

"""
    chebyshev!(x, A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]

Solve Ax = b for symmetric, definite matrices A using Chebyshev iteration.

# Arguments

- `x`: initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side;
- `λmin::Real`: lower bound for the real eigenvalues
- `λmax::Real`: upper bound for the real eigenvalues

## Keywords

- `initially_zero::Bool = false`: if `true` assumes that `iszero(x)` so that one 
  matrix-vector product can be saved when computing the initial 
  residual vector;
- `tol`: tolerance for stopping condition `|r_k| / |r_0| ≤ tol`.
- `maxiter::Int = size(A, 2)`: maximum number of inner iterations of GMRES;
- `Pl = Identity()`: left preconditioner;
- `log::Bool = false`: keep track of the residual norm in each iteration;
- `verbose::Bool = false`: print convergence information during the iterations.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;
- `history`: convergence history.
"""
function chebyshev!(x, A, b, λmin::Real, λmax::Real;
    Pl = Identity(),
    tol::Real=sqrt(eps(real(eltype(b)))),
    maxiter::Int=size(A, 2),
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
