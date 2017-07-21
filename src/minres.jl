export minres_iterable, minres

import Base.LinAlg: BLAS.axpy!, givensAlgorithm
import Base: start, next, done

"""
MINRES is full GMRES for Hermetian matrices, finding
xₖ := x₀ + Vₖyₖ, where Vₖ is an orthonormal basis for
the Krylov subspace K(A, b - Ax₀). Since
|rₖ| = |b - Axₖ| = |r₀ - Vₖ₊₁Hₖyₖ| = | |r₀|e₁ - Hₖyₖ|, we solve
Hₖyₖ = |r₀|e₁ in least-square sense. As the Hessenberg matrix
is tridiagonal, its QR decomp Hₖ = QₖRₖ has Rₖ with only 3 diagonals,
easily obtained with Givens rotations. The least-squares problem is
solved as Hₖ'Hₖyₖ = Hₖ'|r₀|e₁ => yₖ = inv(Rₖ) Qₖ'|r₀|e₁, so
xₖ := x₀ + [Vₖ inv(Rₖ)] [Qₖ'|r₀|e₁]. Now the main difference with GMRES
is the placement of the brackets. MINRES computes Wₖ = Vₖ inv(Rₖ) via 3-term
recurrences using only the last column of R, and computes last 
two terms of Qₖ'|r₀|e₁ as well.
The active column of the Hessenberg matrix H is updated in place to form the
active column of R. 
"""
type MINRESIterable{matT, vecT, realVecT, realT}
    A::matT
    x::vecT

    # Krylov basis vectors
    v_prev::vecT
    v_curr::vecT
    v_next::vecT

    # W = R * inv(V) is computed using 3-term recurrence
    w_prev::vecT
    w_curr::vecT
    w_next::vecT

    # Vector of size 4, holding the active column of the Hessenberg matrix
    # rhs is just two active values of the right-hand side.
    H::realVecT
    rhs::realVecT

    # Some Givens rotations
    c_prev::realT
    s_prev::realT
    c_curr::realT
    s_curr::realT

    # The normalization constant of v_prev
    prev_norm::realT

    # Bookkeeping
    mv_products::Int
    maxiter::Int
    tolerance::realT
    resnorm::realT
end

minres_iterable(A, b; kwargs...) = minres_iterable!(zerox(A, b), A, b; initially_zero = true, kwargs...)

function minres_iterable!(x, A, b; initially_zero = false, tol = sqrt(eps(real(eltype(b)))), maxiter = size(A, 1))
    T = real(eltype(b))

    v_prev = similar(b)
    v_curr = copy(b)
    v_next = similar(b)
    w_prev = similar(b)
    w_curr = similar(b)
    w_next = similar(b)

    mv_products = 0

    # For nonzero x's, we must do an MV for the initial residual vec
    if !initially_zero
        # Use v_next to store Ax; v_next will soon be overwritten.
        A_mul_B!(v_next, A, x)
        axpy!(-one(T), v_next, v_curr)
        mv_products = 1
    end

    # Last active column of the Hessenberg matrix
    H = zeros(T, 4)

    resnorm = norm(v_curr)
    reltol = resnorm * tol

    # Last two entries of the right-hand side
    rhs = [resnorm; zero(T)]

    # Normalize the first Krylov basis vector
    scale!(v_curr, inv(resnorm))

    # The normalization constant of v_prev (initially zero)
    prev_norm = zero(T)

    # Givens rotations
    c_prev, s_prev = one(T), zero(T)
    c_curr, s_curr = one(T), zero(T)

    MINRESIterable(
        A, x,
        v_prev, v_curr, v_next,
        w_prev, w_curr, w_next,
        H, rhs,
        c_prev, s_prev, c_curr, s_curr, prev_norm,
        mv_products, maxiter, reltol, resnorm
    )
end

converged(m::MINRESIterable) = m.resnorm ≤ m.tolerance

start(::MINRESIterable) = 1

done(m::MINRESIterable, iteration::Int) = iteration > m.maxiter || converged(m)

function next(m::MINRESIterable, iteration::Int)
    # v_next = A * v_curr - prev_norm * v_prev
    A_mul_B!(m.v_next, m.A, m.v_curr)

    iteration > 1 && axpy!(-m.prev_norm, m.v_prev, m.v_next)
    
    # Orthogonalize w.r.t. v_curr. The tridiagonal Hessenberg
    # matrix is real, so let's enforce real arithmetic
    m.H[3] = real(dot(m.v_curr, m.v_next))
    axpy!(-m.H[3], m.v_curr, m.v_next)

    # Normalize
    m.H[4] = m.prev_norm = norm(m.v_next)
    scale!(m.v_next, inv(m.prev_norm))

    # Rotation on H[1] and H[2]. Note that H[1] = 0 initially
    if iteration > 2
        m.H[1] = m.s_prev * m.H[2]
        m.H[2] = m.c_prev * m.H[2]
    end

    # Rotation on H[2] and H[3]
    if iteration > 1
        tmp = -m.s_curr * m.H[2] + m.c_curr * m.H[3]
        m.H[2] = m.c_curr * m.H[2] + m.s_curr * m.H[3]
        m.H[3] = tmp
    end

    # Next rotation
    c, s, m.H[3] = givensAlgorithm(m.H[3], m.H[4])

    # Apply as well to the right-hand side
    m.rhs[2] = -s * m.rhs[1]
    m.rhs[1] = c * m.rhs[1]

    # Update W = V * inv(R). Two axpy's can maybe be one MV.
    copy!(m.w_next, m.v_curr)
    iteration > 1 && axpy!(-m.H[2], m.w_curr, m.w_next)
    iteration > 2 && axpy!(-m.H[1], m.w_prev, m.w_next)
    scale!(m.w_next, inv(m.H[3]))

    # Update solution x
    axpy!(m.rhs[1], m.w_next, m.x)

    # Move on: next -> curr, curr -> prev
    m.v_prev, m.v_curr, m.v_next = m.v_curr, m.v_next, m.v_prev
    m.w_prev, m.w_curr, m.w_next = m.w_curr, m.w_next, m.w_prev
    m.c_prev, m.s_prev, m.c_curr, m.s_curr = m.c_curr, m.s_curr, c, s
    m.rhs[1] = m.rhs[2]

    # Due to symmetry of the tri-diagonal matrix
    m.H[2] = m.prev_norm

    # The approximate residual is cheaply available
    m.resnorm = abs(m.rhs[2])

    m.resnorm, iteration + 1
end

function minres!(x, A, b; 
    verbose::Bool = false,
    log::Bool = false,
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Int = min(30, size(A, 1)),
    kwargs...
)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && reserve!(history, :resnorm, maxiter)
    
    iterable = minres_iterable!(x, A, b; tol = tol, maxiter = maxiter, kwargs...)
    
    if log
        history.mvps = iterable.mv_products
    end

    for (iteration, resnorm) = enumerate(iterable)
        if log
            nextiter!(history, mvps = 1)
            push!(history, :resnorm, resnorm)
        end
        verbose && @printf("%3d\t%1.2e\n", iteration, resnorm)
    end
    
    verbose && println()
    log && setconv(history, converged(iterable))
    log && shrink!(history)
    
    log ? (iterable.x, history) : iterable.x
end

minres(A, b; kwargs...) = minres!(zerox(A, b), A, b; initially_zero = true, kwargs...)

#################
# Documentation #
#################

let
doc_call = "minres(A, b)"
doc!_call = "minres!(x, A, b)"

doc_msg = "Using initial guess zeros(b)."
doc!_msg = "Overwrites `x`."

doc_arg = ""
doc!_arg = """`x`: initial guess, overwrite final approximation."""

doc_version = (doc_call, doc_msg, doc_arg)
doc!_version = (doc!_call, doc!_msg, doc!_arg)

docstring = String[]

#Build docs
for (call, msg, arg) in (doc_version, doc!_version) #Start
    push!(docstring, 
"""
$call

Solve A*x = b for Hermetian matrices A using MINRES. The method is mathematically 
equivalent to unrestarted GMRES, but exploits symmetry of A, resulting in short
recurrences requiring only 6 vectors of storage. MINRES might be slightly less
stable than full GMRES.

$msg

# Arguments

$arg

`A`: linear operator.

`b`: right hand side (vector).

## Keywords

`tol::Real = sqrt(eps(real(eltype(b))))`: tolerance for stopping condition 
`|r_k| / |r_0| ≤ tol`. Note that the residual is computed only approximately.

`maxiter::Int = min(30, size(A, 1))`: maximum number of iterations.

`verbose::Bool = false` output during the iterations

`log::Bool = false` enables logging, see **Output**.

# Output

**if `log` is `false`**

`x`: approximated solution.

**if `log` is `true`**

`x`: approximated solution.

`ch`: convergence history.
"""
    )
end

@doc docstring[1] -> minres
@doc docstring[2] -> minres!
end