import Base: iterate
using Printf
export cocg, cocg!, COCGIterable, PCOCGIterable, cocg_iterator!, COCGStateVariables

mutable struct COCGIterable{matT, solT, vecT, numT <: Real, paramT <: Number}
    A::matT
    x::solT
    r::vecT
    c::vecT
    u::vecT
    tol::numT
    residual::numT
    rr_prev::paramT
    maxiter::Int
    mv_products::Int
end

mutable struct PCOCGIterable{precT, matT, solT, vecT, numT <: Real, paramT <: Number}
    Pl::precT
    A::matT
    x::solT
    r::vecT
    c::vecT
    u::vecT
    tol::numT
    residual::numT
    rc_prev::paramT
    maxiter::Int
    mv_products::Int
end

@inline converged(it::Union{COCGIterable, PCOCGIterable}) = it.residual ≤ it.tol

@inline start(it::Union{COCGIterable, PCOCGIterable}) = 0

@inline done(it::Union{COCGIterable, PCOCGIterable}, iteration::Int) = iteration ≥ it.maxiter || converged(it)


#################
# Ordinary COCG #
#################

function iterate(it::COCGIterable, iteration::Int=start(it))
    # Check for termination first
    if done(it, iteration)
        return nothing
    end

    # u := r + βu (almost an axpy)
    rr = sum(rₖ^2 for rₖ in it.r)  # rᵀ * r
    β = rr / it.rr_prev
    it.u .= it.r .+ β .* it.u

    # c = A * u
    mul!(it.c, it.A, it.u)
    uc = sum(uₖ*cₖ for (uₖ,cₖ) in zip(it.u,it.c))  # uᵀ * c
    α = rr / uc

    # Improve solution and residual
    it.rr_prev = rr
    it.x .+= α .* it.u
    it.r .-= α .* it.c

    it.residual = norm(it.r)

    # Return the residual at item and iteration number as state
    it.residual, iteration + 1
end

#######################
# Preconditioned COCG #
#######################

function iterate(it::PCOCGIterable, iteration::Int=start(it))
    # Check for termination first
    if done(it, iteration)
        return nothing
    end

    # Apply left preconditioner: c = Pl⁻¹ r
    ldiv!(it.c, it.Pl, it.r)

    # u := c + βu (almost an axpy)
    rc = sum(rₖ*cₖ for (rₖ,cₖ) in zip(it.r,it.c))  # rᵀ * c
    β = rc / it.rc_prev
    it.u .= it.c .+ β .* it.u

    # c = A * u
    mul!(it.c, it.A, it.u)
    uc = sum(uₖ*cₖ for (uₖ,cₖ) in zip(it.u,it.c))  # uᵀ * c
    α = rc / uc

    # Improve solution and residual
    it.rc_prev = rc
    it.x .+= α .* it.u
    it.r .-= α .* it.c

    it.residual = norm(it.r)

    # Return the residual at item and iteration number as state
    it.residual, iteration + 1
end

# Utility functions

"""
Intermediate COCG state variables to be used inside cocg and cocg!. `u`, `r` and `c` should be of the same type as the solution of `cocg` or `cocg!`.
```
struct COCGStateVariables{T,Tx<:AbstractArray{T}}
    u::Tx
    r::Tx
    c::Tx
end
```
"""
struct COCGStateVariables{T,Tx<:AbstractArray{T}}
    u::Tx
    r::Tx
    c::Tx
end

function cocg_iterator!(x, A, b, Pl = Identity();
                        abstol::Real = zero(real(eltype(b))),
                        reltol::Real = sqrt(eps(real(eltype(b)))),
                        maxiter::Int = size(A, 2),
                        statevars::COCGStateVariables = COCGStateVariables(zero(x), similar(x), similar(x)),
                        initially_zero::Bool = false)
    u = statevars.u
    r = statevars.r
    c = statevars.c
    u .= zero(eltype(x))
    copyto!(r, b)

    # Compute r with an MV-product or not.
    if initially_zero
        mv_products = 0
    else
        mv_products = 1
        mul!(c, A, x)
        r .-= c
    end
    residual = norm(r)
    tolerance = max(reltol * residual, abstol)

    # Return the iterable
    if isa(Pl, Identity)
        return COCGIterable(A, x, r, c, u,
            tolerance, residual, one(eltype(r)),
            maxiter, mv_products
        )
    else
        return PCOCGIterable(Pl, A, x, r, c, u,
            tolerance, residual, one(eltype(r)),
            maxiter, mv_products
        )
    end
end

"""
    cocg(A, b; kwargs...) -> x, [history]

Same as [`cocg!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
cocg(A, b; kwargs...) = cocg!(zerox(A, b), A, b; initially_zero = true, kwargs...)

"""
    cocg!(x, A, b; kwargs...) -> x, [history]

# Arguments

- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `statevars::COCGStateVariables`: Has 3 arrays similar to `x` to hold intermediate results;
- `initially_zero::Bool`: If `true` assumes that `iszero(x)` so that one
  matrix-vector product can be saved when computing the initial
  residual vector;
- `Pl = Identity()`: left preconditioner of the method. Should be symmetric,
  positive-definite like `A`;
- `abstol::Real = zero(real(eltype(b)))`,
  `reltol::Real = sqrt(eps(real(eltype(b))))`: absolute and relative
  tolerance for the stopping condition
  `|r_k| / |r_0| ≤ max(reltol * resnorm, abstol)`, where `r_k = A * x_k - b`
  is the residual in the `k`th iteration;
- `maxiter::Int = size(A,2)`: maximum number of iterations;
- `verbose::Bool = false`: print method information;
- `log::Bool = false`: keep track of the residual norm in each iteration.

# Output

**if `log` is `false`**

- `x`: approximated solution.

**if `log` is `true`**

- `x`: approximated solution.
- `ch`: convergence history.

**ConvergenceHistory keys**

- `:tol` => `::Real`: stopping tolerance.
- `:resnom` => `::Vector`: residual norm at each iteration.
"""
function cocg!(x, A, b;
               abstol::Real = zero(real(eltype(b))),
               reltol::Real = sqrt(eps(real(eltype(b)))),
               maxiter::Int = size(A, 2),
               log::Bool = false,
               statevars::COCGStateVariables = COCGStateVariables(zero(x), similar(x), similar(x)),
               verbose::Bool = false,
               Pl = Identity(),
               kwargs...)
    history = ConvergenceHistory(partial = !log)
    history[:abstol] = abstol
    history[:reltol] = reltol
    log && reserve!(history, :resnorm, maxiter + 1)

    # Actually perform COCG
    iterable = cocg_iterator!(x, A, b, Pl; abstol = abstol, reltol = reltol, maxiter = maxiter,
                              statevars = statevars, kwargs...)
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

    log ? (iterable.x, history) : iterable.x
end
