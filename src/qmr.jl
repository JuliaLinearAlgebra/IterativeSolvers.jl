import Base: iterate

export qmr, qmr!

mutable struct LanczosDecomp{T, rT, matT, mat_tT, vecT}
    A::matT
    At::mat_tT

    v_prev::vecT # Orthonormal basis vectors for A
    v_curr::vecT # Orthonormal basis vectors for A
    v_next::vecT # Orthonormal basis vectors for A

    w_prev::vecT # Orthonormal basis vectors for A'
    w_curr::vecT # Orthonormal basis vectors for A'
    w_next::vecT # Orthonormal basis vectors for A'

    α::T
    β_prev::T
    β_curr::T
    δ::T
    resnorm::rT
end
function LanczosDecomp(x, A::matT, b; initially_zero = false) where {matT}
    T = eltype(A)

    # we choose:
    # v_0 = 0, v_1 = b - Ax
    # w_0 = 0, w_1 = v_1
    # v_2 and w_2 are uninitialized and used as workspace

    v_prev = zero(x)
    v_curr = copy(b)
    v_next = similar(x)

    if !initially_zero
        # r0 = A*x - b
        mul!(v_next, A, x)
        axpy!(-one(T), v_next, v_curr)
    end
    resnorm = norm(v_curr)
    rmul!(v_curr, inv(resnorm))

    w_prev = zero(x)
    w_curr = copy(v_curr)
    w_next = similar(x)

    α = zero(T)
    β_prev = zero(T)
    β_curr = zero(T)
    δ = zero(T)

    LanczosDecomp(
        A, adjoint(A),
        v_prev, v_curr, v_next,
        w_prev, w_curr, w_next,
        α, β_prev, β_curr, δ,
        resnorm)
end

struct LookAheadLanczosDecomp{T, matT, mat_tT, vecT}
    A::matT
    At::mat_tT
    Ax::vecT

    v_prev::matT # Orthonormal basis vectors for A
    v_curr::matT # Orthonormal basis vectors for A
    v_next::matT # Orthonormal basis vectors for A

    w_prev::matT # Orthonormal basis vectors for A'
    w_curr::matT # Orthonormal basis vectors for A'
    w_next::matT # Orthonormal basis vectors for A'

    δ_curr::matT # adjoint(w)* v
    δ_prev::matT

    WtAv_curr::matT
    WtAv_prev::matT

    H::Matrix{T} # Hessenberg Matrix, computed on CPU
end
function LookAheadLanczosDecomp(A::matT, b) where {matT}
    # TODO: implement
    # v̂_(n+1) = Av_n - V_l δ_l \ W_l' A v_n - V_l-1 δ_(l-1) W_(l-1)' A v_n
    # ŵ_(n+1) = A'w_n - W_l δ_l' \ V_l' A' w_n - W_l-1 δ_(l-1)' V_(l-1)' A' v_n
end

start(::LanczosDecomp) = 1
done(l::LanczosDecomp, iteration::Int) = false
function iterate(l::LanczosDecomp, iteration::Int=start(l))
    # Following Algorithm 7.1 of Saad
    if done(l, iteration) return nothing end

    mul!(l.v_next, l.A, l.v_curr)

    l.α = dot(l.v_next, l.w_curr) # v_next currently contains A*v_curr
    axpy!(-conj(l.α), l.v_curr, l.v_next)
    if iteration > 1 # β_1 = 0
        axpy!(-conj(l.β_curr), l.v_prev, l.v_next)
    end

    mul!(l.w_next, l.At, l.w_curr)
    axpy!(-l.α, l.w_curr, l.w_next)
    if iteration > 1 # δ_1 = 0
        axpy!(-l.δ, l.w_prev, l.w_next)
    end

    vw = dot(l.v_next, l.w_next)
    l.δ = sqrt(abs(vw))
    if iszero(l.δ) return nothing end

    l.β_prev = l.β_curr
    l.β_curr = vw / l.δ

    rmul!(l.v_next, inv(l.δ))
    rmul!(l.w_next, inv(l.β_curr))

    l.w_prev .= l.w_curr
    l.w_curr .= l.w_next
    l.v_prev .= l.v_curr
    l.v_curr .= l.v_next

    return nothing, iteration + 1
end


mutable struct QMRIterable{T, xT, rT, preclT, precrT, lanczosT}
    x::xT

    Pl::preclT
    Pr::precrT
    lanczos::lanczosT
    resnorm::rT
    reltol::rT
    maxiter::Int

    g::Vector{T} # right-hand size following Hessenberg multiplication
    H::Vector{T}  # Hessenberg Matrix, computed on CPU

    c_prev::T
    c_curr::T
    s_prev::T
    s_curr::T

    p_prev::xT
    p_curr::xT
end

function QMRIterable(x, A, b;
    Pl = Identity(),
    Pr = Identity(),
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Int = size(A, 2),
    initially_zero::Bool = false,
    lookahead::Bool = false
    )
    T = eltype(x)

    # Approximate solution
    if lookahead
        lanczos = LookAheadLanczosDecomp(x, A, b, initially_zero = initially_zero)
    else
        lanczos = LanczosDecomp(x, A, b, initially_zero = initially_zero)
    end
    # A' = Pl \ A / Pr
    # b' = Pl \ (b - Ax_0)
    # A'y = b'
    # y  = Pr * (x - x_0)
    # x_n = x_0 + Pr \ y
    # r_n = Pl r'_n

    resnorm = lanczos.resnorm
    g = [resnorm; zero(T)]
    H = zeros(T, 4)

    # Givens rotations
    c_prev, s_prev = one(T), zero(T)
    c_curr, s_curr = one(T), zero(T)

    # Projection operators
    p_prev = zero(x)
    p_curr = zero(x)

    reltol = tol * lanczos.resnorm

    QMRIterable(x,
        Pl, Pr,
        lanczos, resnorm, reltol,
        maxiter,
        g, H,
        c_prev, c_curr, s_prev, s_curr,
        p_prev, p_curr
    )
end

converged(q::QMRIterable) = q.resnorm ≤ q.reltol
start(::QMRIterable) = 1
done(q::QMRIterable, iteration::Int) = iteration > q.maxiter || converged(q)

function iterate(q::QMRIterable, iteration::Int=start(q))
    if done(q, iteration) return nothing end

    iterate(q.lanczos, iteration)

    # for classical Lanczos algorithm, only need 4 elements of last column
    # of Hessenberg matrix
    # H[1] is h_m,(m-2), H[2] is h_m,(m-1), H[3] is h_m,m, and H[4] is h_m,(m+1)
    q.H[2] = conj(q.lanczos.β_prev) # β_m
    q.H[3] = conj(q.lanczos.α)      # α_m
    q.H[4] = q.lanczos.δ            # δ_(m+1)

    # Rotation on H[1] and H[2]. Note that H[1] = 0 initially
    if iteration > 2
        q.H[1] = q.s_prev * q.H[2]
        q.H[2] = q.c_prev * q.H[2]
    end

    # Rotation on H[2] and H[3]
    if iteration > 1
        tmp    = -conj(q.s_curr) * q.H[2] + q.c_curr * q.H[3]
        q.H[2] =       q.c_curr  * q.H[2] + q.s_curr * q.H[3]
        q.H[3] = tmp
    end

    a, b = q.H[3], q.H[4]
    # Compute coefficients for Ω_m, and apply it
    c, s, q.H[3] = givensAlgorithm(q.H[3], q.H[4])

    # Apply Ω_m to rhs
    q.g[2] = -conj(s) * q.g[1]
    q.g[1] =       c  * q.g[1]

    # H is now properly factorized
    # compute p_m
    # we re-use q.lanczos.v_next as workspace for p_m
    q.lanczos.v_next .= q.lanczos.v_prev # we need v_m, not v_m+1
    iteration > 1 && axpy!(-q.H[2], q.p_curr, q.lanczos.v_next)
    iteration > 2 && axpy!(-q.H[1], q.p_prev, q.lanczos.v_next)
    rmul!(q.lanczos.v_next, inv(q.H[3]))

    # iterate x_n = x_n-1 + γ_m p_m
    axpy!(q.g[1], q.lanczos.v_next, q.x)

    # transfer coefficients of Givens rotations for next iteration
    q.c_prev, q.s_prev, q.c_curr, q.s_curr = q.c_curr, q.s_curr, c, s
    q.p_prev .= q.p_curr
    q.p_curr .= q.lanczos.v_next
    q.g[1] = q.g[2]

    # approximate residual
    # See Proposition 7.3 of Saad
    q.resnorm = abs(q.g[2])

    q.resnorm, iteration + 1
end

"""
    qmr(A, b; kwargs...) -> x, [history]

Same as [`qmr!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
qmr(A, b; kwargs...) = qmr!(zerox(A, b), A, b; initially_zero = true, kwargs...)

"""
    qmr!(x, A, b; kwargs...) -> x, [history]

Solves the problem ``Ax = b`` with the Quasi-Minimal Residual (QMR) method.

# Arguments
- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `initally_zero::Bool`: If `true` assumes that `iszero(x)` so that one
  matrix-vector product can be saved when computing the initial residual
  vector;
- `maxiter::Int = size(A, 2)`: maximum number of iterations;
- `tol`: relative tolerance;
- `Pl`: left preconditioner;
- `Pr`: right preconditioner;
- `lookahead::Bool`: use look-ahead Lanczos process
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iteration.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;

- `history`: convergence history.
"""
function qmr!(x, A, b;
    Pl = Identity(),
    Pr = Identity(),
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Int = size(A, 2),
    lookahead::Bool = false,
    log::Bool = false,
    initially_zero::Bool = false,
    verbose::Bool = false
)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && reserve!(history, :resnorm, maxiter)

    iterable = QMRIterable(x, A, b; Pl = Pl, Pr = Pr, tol = tol, maxiter = maxiter, initially_zero = initially_zero)

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
