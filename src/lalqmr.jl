import Base: iterate

export lalqmr, lalqmr!

mutable struct LALQMRIterable{T, xT, MT, rT, lanczosT}
    x::xT

    lanczos::lanczosT
    resnorm::rT
    tol::rT
    maxiter::Int

    t::Vector{T}
    R::LimitedMemoryUpperTriangular{T, Matrix{T}}

    G::Vector{Givens{T}}

    d::LimitedMemoryMatrix{T, MT}
    dlast::xT
end

function lalqmr_iterable!(
    x, A, b;
    abstol::Real = zero(real(eltype(b))),
    reltol::Real = sqrt(eps(real(eltype(b)))),
    maxiter::Int = size(A, 2),
    initially_zero::Bool = false,
    max_memory::Int=4,
    kwargs...
)
    T = eltype(x)
    r = copy(b)
    if !initially_zero
        r -= A*x
    end
    resnorm = norm(r)
    v = r/resnorm
    w = copy(v)

    lanczos = LookAheadLanczosDecomp(
        A, v, w;
        vw_normalized=true,
        max_memory=max_memory,
        kwargs...
    )

    R = LimitedMemoryUpperTriangular{T, Matrix{T}}(max_memory)
    # Givens rotations
    G = Vector{Givens{T}}()

    # Projection operators
    d = LimitedMemoryMatrix(similar(x, size(v, 1), 0), max_memory)

    t = Vector{T}(undef, 1)
    t[1] = resnorm

    tolerance = max(reltol * resnorm, abstol)
    dlast = similar(x)

    return LALQMRIterable(
        x,
        lanczos, resnorm, tolerance,
        maxiter,
        t, R,
        G, d, dlast
    )
end

converged(q::LALQMRIterable) = q.resnorm ≤ q.tol
start(::LALQMRIterable) = 1
done(q::LALQMRIterable, iteration::Int) = iteration > q.maxiter || converged(q)

function iterate(q::LALQMRIterable, n::Int=start(q))
    # Check for termination first
    if done(q, n)
        return nothing
    end

    iterate(q.lanczos, n)

    # Eq. 6.2, update QR factorization of L
    Llastcol = q.lanczos.L[:, end]
    for g in q.G
        Llastcol = g*Llastcol
    end
    gend, r = givens(Llastcol, n, n+1)
    push!(q.G, gend)
    Llastcol[end-1] = r # Llastcol[end] = 0, but we don't need it
    _grow_hcat!(q.R, Llastcol[1:end-1])

    # Eq. 6.2, update t
    push!(q.t, 0)
    q.t .= gend * q.t
  
    # Eq. 6.3, calculate projection
    # Dn = [Vn Un^-1 Rn^-1]
    # => Dn Rn Un = Vn
    RU = q.R*q.lanczos.U
    copyto!(q.dlast, view(q.lanczos.V, :, n))
    for i in 1:size(RU, 1)-1
        if RU[i, end] != 0
            axpy!(-RU[i, end], view(q.d, :, i), q.dlast)
        end
    end
    q.dlast ./= RU[end, end]
    hcat!(q.d, q.dlast)

    # iterate x_n = x_n-1 + d_n τ_n
    axpy!(q.t[end-1], q.dlast, q.x)

    # Eq. 4.12, Freund 1990
    q.resnorm = q.resnorm * abs(gend.s) * sqrt(n+1)/sqrt(n)

    return q.resnorm, n + 1
end

"""
    lalqmr(A, b; kwargs...) -> x, [history]

Same as [`lalqmr!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
lalqmr(A, b; kwargs...) = lalqmr!(zerox(A, b), A, b; initially_zero = true, kwargs...)

"""
    lalqmr!(x, A, b; kwargs...) -> x, [history]

Solves the problem ``Ax = b`` with the Quasi-Minimal Residual (QMR) method with Look-ahead. See [`LookAheadLanczosDecomp`](@ref).

# Arguments
- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `initally_zero::Bool`: If `true` assumes that `iszero(x)` so that one
  matrix-vector product can be saved when computing the initial residual
  vector;
- `maxiter::Int = size(A, 2)`: maximum number of iterations;
- `abstol::Real = zero(real(eltype(b)))`,
  `reltol::Real = sqrt(eps(real(eltype(b))))`: absolute and relative
  tolerance for the stopping condition
  `|r_k| ≤ max(reltol * |r_0|, abstol)`, where `r_k = A * x_k - b`
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iteration.
- `max_block_size`: maximum block size during look-ahead process
- `max_memory`: maximum allowed memory for look-ahead process

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;

- `history`: convergence history.

[^Freund1990]:
    Freund, W. R., & Nachtigal, N. M. (1990). QMR : for a Quasi-Minimal
    Residual Linear Method Systems. (December).
"""
function lalqmr!(
    x, A, b;
    abstol::Real = zero(real(eltype(b))),
    reltol::Real = sqrt(eps(real(eltype(b)))),
    maxiter::Int = size(A, 2),
    log::Bool = false,
    initially_zero::Bool = false,
    verbose::Bool = false,
    kwargs...
)
    history = ConvergenceHistory(partial = !log)
    history[:abstol] = abstol
    history[:reltol] = reltol
    log && reserve!(history, :resnorm, maxiter)

    iterable = lalqmr_iterable!(
        x, A, b;
        abstol=abstol,
        reltol=reltol,
        maxiter=maxiter,
        initially_zero=initially_zero,
        log=log,
        verbose=verbose,
        kwargs...
    )

    verbose && @printf("=== qmr ===\n%4s\t%7s\n","iter","resnorm")

    for (iteration, residual) = enumerate(iterable)
        if log
            nextiter!(history)
            push!(history, :resnorm, residual)
        end

        verbose && @printf("%3d\t%1.2e\n", iteration, residual)
    end

    verbose && println()
    log && setconv(history, converged(iterable))
    log && shrink!(history)

    log ? (x, history) : x
end