export bicgstabl, bicgstabl!, bicgstabl_iterator, bicgstabl_iterator!, BiCGStabIterable

import Base: start, next, done

mutable struct BiCGStabIterable{precT, matT, solT, vecT <: AbstractVector, smallMatT <: AbstractMatrix, realT <: Real, scalarT <: Number}
    A::matT
    l::Int

    x::solT
    r_shadow::vecT
    rs::smallMatT
    us::smallMatT

    max_mv_products::Int
    mv_products::Int
    reltol::realT
    residual::realT

    Pl::precT

    γ::vecT
    ω::scalarT
    σ::scalarT
    M::smallMatT
end

bicgstabl_iterator(A, b, l; kwargs...) = bicgstabl_iterator!(zerox(A, b), A, b, l; initial_zero = true, kwargs...)

function bicgstabl_iterator!(x, A, b, l::Int = 2;
    Pl = Identity(),
    max_mv_products = min(30, size(A, 1)),
    initial_zero = false,
    tol = sqrt(eps(real(eltype(b))))
)
    T = eltype(x)
    n = size(A, 1)
    mv_products = 0

    # Large vectors.
    r_shadow = rand(T, n)
    rs = Matrix{T}(n, l + 1)
    us = zeros(T, n, l + 1)

    residual = view(rs, :, 1)
    
    # Compute the initial residual rs[:, 1] = b - A * x
    # Avoid computing A * 0.
    if initial_zero
        copy!(residual, b)
    else
        A_mul_B!(residual, A, x)
        @blas! residual -= one(T) * b
        @blas! residual *= -one(T)
        mv_products += 1
    end

    # Apply the left preconditioner
    A_ldiv_B!(Pl, residual)

    γ = zeros(T, l)
    ω = σ = one(T)

    nrm = norm(residual)

    # For the least-squares problem
    M = zeros(T, l + 1, l + 1)

    # Stopping condition based on relative tolerance.
    reltol = nrm * tol

    BiCGStabIterable(A, l, x, r_shadow, rs, us,
        max_mv_products, mv_products, reltol, nrm,
        Pl,
        γ, ω, σ, M
    )
end

@inline converged(it::BiCGStabIterable) = it.residual ≤ it.reltol
@inline start(::BiCGStabIterable) = 0
@inline done(it::BiCGStabIterable, iteration::Int) = it.mv_products ≥ it.max_mv_products || converged(it)

function next(it::BiCGStabIterable, iteration::Int)
    T = eltype(it.x)
    L = 2 : it.l + 1

    it.σ = -it.ω * it.σ
    
    ## BiCG part
    for j = 1 : it.l
        ρ = dot(it.r_shadow, view(it.rs, :, j))
        β = ρ / it.σ
        
        # us[:, 1 : j] .= rs[:, 1 : j] - β * us[:, 1 : j]
        for i = 1 : j
            @blas! view(it.us, :, i) *= -β
            @blas! view(it.us, :, i) += one(T) * view(it.rs, :, i)
        end

        # us[:, j + 1] = Pl \ (A * us[:, j])
        next_u = view(it.us, :, j + 1)
        A_mul_B!(next_u, it.A, view(it.us, :, j))
        A_ldiv_B!(it.Pl, next_u)

        it.σ = dot(it.r_shadow, next_u)
        α = ρ / it.σ

        # rs[:, 1 : j] .= rs[:, 1 : j] - α * us[:, 2 : j + 1]
        for i = 1 : j
            @blas! view(it.rs, :, i) -= α * view(it.us, :, i + 1)
        end
        
        # rs[:, j + 1] = Pl \ (A * rs[:, j])
        next_r = view(it.rs, :, j + 1)
        A_mul_B!(next_r, it.A , view(it.rs, :, j))
        A_ldiv_B!(it.Pl, next_r)
        
        # x = x + α * us[:, 1]
        @blas! it.x += α * view(it.us, :, 1)
    end

    # Bookkeeping
    it.mv_products += 2 * it.l

    ## MR part
    
    # M = rs' * rs
    Ac_mul_B!(it.M, it.rs, it.rs)

    # γ = M[L, L] \ M[L, 1] 
    F = lufact!(view(it.M, L, L))
    A_ldiv_B!(it.γ, F, view(it.M, L, 1))

    # This could even be BLAS 3 when combined.
    BLAS.gemv!('N', -one(T), view(it.us, :, L), it.γ, one(T), view(it.us, :, 1))
    BLAS.gemv!('N', one(T), view(it.rs, :, 1 : it.l), it.γ, one(T), it.x)
    BLAS.gemv!('N', -one(T), view(it.rs, :, L), it.γ, one(T), view(it.rs, :, 1))

    it.ω = it.γ[it.l]
    it.residual = norm(view(it.rs, :, 1))

    it.residual, iteration + 1
end

# Classical API

"""
    bicgstabl(A, b, l; kwargs...) -> x, [history]

Same as [`bicgstabl!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
bicgstabl(A, b, l = 2; kwargs...) = bicgstabl!(zerox(A, b), A, b, l; initial_zero = true, kwargs...)

"""
    bicgstabl!(x, A, b, l; kwargs...) -> x, [history]

# Arguments

- `A`: linear operator;
- `b`: right hand side (vector);
- `l::Int = 2`: Number of GMRES steps.

## Keywords

- `max_mv_products::Int = min(30, size(A, 1))`: maximum number of matrix vector products.
For BiCGStab(l) this is a less dubious term than "number of iterations";
- `Pl = Identity()`: left preconditioner of the method;
- `tol::Real = sqrt(eps(real(eltype(b))))`: tolerance for stopping condition `|r_k| / |r_0| ≤ tol`. 
   Note that (1) the true residual norm is never computed during the iterations, 
   only an approximation; and (2) if a preconditioner is given, the stopping condition is based on the 
   *preconditioned residual*.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;
- `history`: convergence history.
"""
function bicgstabl!(x, A, b, l = 2;
    tol = sqrt(eps(real(eltype(b)))),
    max_mv_products::Int = min(20, size(A, 1)),
    log::Bool = false,
    verbose::Bool = false,
    Pl = Identity(),
    kwargs...
)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol

    # This doesn't yet make sense: the number of iters is smaller.
    log && reserve!(history, :resnorm, max_mv_products)
    
    # Actually perform CG
    iterable = bicgstabl_iterator!(x, A, b, l; Pl = Pl, tol = tol, max_mv_products = max_mv_products, kwargs...)
    
    if log
        history.mvps = iterable.mv_products
    end

    for (iteration, item) = enumerate(iterable)
        if log
            nextiter!(history)
            history.mvps = iterable.mv_products
            push!(history, :resnorm, iterable.residual)
        end
        verbose && @printf("%3d\t%1.2e\n", iteration, iterable.residual)
    end
    
    verbose && println()
    log && setconv(history, converged(iterable))
    log && shrink!(history)
    
    log ? (iterable.x, history) : iterable.x
end
