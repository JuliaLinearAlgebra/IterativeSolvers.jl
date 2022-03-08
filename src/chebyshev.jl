import Base: iterate

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
    tol::realT
    maxiter::Int
    mv_products::Int
end

@inline converged(it::ChebyshevIterable) = it.resnorm ≤ it.tol
@inline start(::ChebyshevIterable) = 0
@inline done(it::ChebyshevIterable, iteration::Int) = iteration ≥ it.maxiter || converged(it)

function iterate(cheb::ChebyshevIterable, iteration::Int=start(cheb))
    # Check for termination first
    if done(cheb, iteration)
        return nothing
    end

    T = eltype(cheb.x)

    ldiv!(cheb.c, cheb.Pl, cheb.r)

    if iteration == 1
        cheb.α = T(2) / cheb.λ_avg
        copyto!(cheb.u, cheb.c)
    else
        β = (cheb.λ_diff * cheb.α / 2) ^ 2
        cheb.α = inv(cheb.λ_avg - β)
        cheb.u .= cheb.c .+ β .* cheb.c
    end

    mul!(cheb.c, cheb.A, cheb.u)
    cheb.mv_products += 1

    cheb.x .+= cheb.α .* cheb.u
    cheb.r .-= cheb.α .* cheb.c

    cheb.resnorm = norm(cheb.r)

    cheb.resnorm, iteration + 1
end

function chebyshev_iterable!(x, A, b, λmin::Real, λmax::Real;
                             abstol::Real = zero(real(eltype(b))),
                             reltol::Real = sqrt(eps(real(eltype(b)))),
                             maxiter = size(A, 2),
                             Pl = Identity(),
                             initially_zero = false)

    λ_avg = (λmax + λmin) / 2
    λ_diff = (λmax - λmin) / 2

    T = eltype(x)
    r = similar(x)
    copyto!(r, b)
    u = zero(x)
    c = similar(x)

    # One MV product less
    if initially_zero
        mv_products = 0
    else
        mv_products = 1
        mul!(c, A, x)
        r .-= c
    end
    resnorm = norm(r)
    tolerance = max(reltol * resnorm, abstol)

    ChebyshevIterable(Pl, A, x, r, u, c,
        zero(real(T)),
        λ_avg, λ_diff,
        resnorm, tolerance, maxiter, mv_products
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
- `abstol::Real = zero(real(eltype(b)))`,
  `reltol::Real = sqrt(eps(real(eltype(b))))`: absolute and relative
  tolerance for the stopping condition
  `|r_k| ≤ max(reltol * |r_0|, abstol)`, where `r_k = A * x_k - b`
  is the residual in the `k`th iteration;
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
                    abstol::Real = zero(real(eltype(b))),
                    reltol::Real = sqrt(eps(real(eltype(b)))),
                    Pl = Identity(),
                    maxiter::Int=size(A, 2),
                    log::Bool=false,
                    verbose::Bool=false,
                    initially_zero::Bool=false)
    history = ConvergenceHistory(partial=!log)
    history[:abstol] = abstol
    history[:reltol] = reltol
    reserve!(history, :resnorm, maxiter)

    verbose && @printf("=========== chebyshev ============\n%4s\t%9s\t%9s\n", "iter", "resnorm", "relresn")

    iterable = chebyshev_iterable!(x, A, b, λmin, λmax; abstol=abstol, reltol=reltol,
                                   maxiter=maxiter, Pl=Pl, initially_zero=initially_zero)
    history.mvps = iterable.mv_products
    resnorm0 = iterable.resnorm
    for (iteration, resnorm) = enumerate(iterable)
        nextiter!(history)
        history.mvps = iterable.mv_products
        push!(history, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.4e\t%1.4e\n", iteration, resnorm, resnorm/resnorm0)
    end
    verbose && println()
    setconv(history, converged(iterable))
    log && shrink!(history)

    log ? (x, history) : x
end
