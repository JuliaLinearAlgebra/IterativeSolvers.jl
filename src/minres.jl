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
recurrences, and computes last two terms of Qₖ'|r₀|e₁ as well. That is
enough to update xₖ each iteration.
"""
type MINRESIterable{matT, vecT, realT}
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
    R::vecT
    rhs::vecT

    # Some Givens rotations
    c_prev::realT
    s_prev::realT
    c_curr::realT
    s_curr::realT

    # The normalization constant of v_prev
    prev_norm::realT

    # Bookkeeping
    maxiter::Int
    tolerance::realT
    resnorm::realT
end

minres_iterable(A, b; kwargs...) = minres_iterable!(zeros(b), A, b; initially_zero = true, kwargs...)

function minres_iterable!(x, A, b; initially_zero = false, tol = sqrt(eps(real(eltype(b)))), maxiter = size(A, 1))
    T = real(eltype(b))

    v_prev = similar(b)
    v_curr = copy(b)
    v_next = similar(b)
    w_prev = similar(b)
    w_curr = similar(b)
    w_next = similar(b)

    # For nonzero x's, we must do an MV for the initial residual vec
    if !initially_zero
        # Use v_next to store Ax; v_next will soon be overwritten.
        A_mul_B!(v_next, A, x)
        axpy!(-one(T), v_next, v_curr)
    end

    # Last column of the R matrix: QR = H where H is 
    # the tridiagonal Lanczos matrix.
    R = zeros(T, 4)

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
        R, rhs,
        c_prev, s_prev, c_curr, s_curr, prev_norm,
        maxiter, reltol, resnorm
    )
end

converged(m::MINRESIterable) = m.resnorm ≤ m.tolerance

start(::MINRESIterable) = 1

done(m::MINRESIterable, iteration::Int) = iteration > m.maxiter || converged(m)

function next(m::MINRESIterable, iteration::Int)
    # v_next = A * v_curr - prev_norm * v_prev
    A_mul_B!(m.v_next, m.A, m.v_curr)

    iteration > 1 && axpy!(-m.prev_norm, m.v_prev, m.v_next)
    
    # Orthogonalize w.r.t. v_curr
    m.R[3] = dot(m.v_curr, m.v_next)
    axpy!(-m.R[3], m.v_curr, m.v_next)

    # Normalize
    m.R[4] = m.prev_norm = norm(m.v_next)
    scale!(m.v_next, inv(m.prev_norm))

    # Rotation on R[1] and R[2]. Note that R[1] = 0 initially
    if iteration > 2
        m.R[1] = m.s_prev * m.R[2]
        m.R[2] = m.c_prev * m.R[2]
    end

    # Rotation on R[2] and R[3]
    if iteration > 1
        tmp = -m.s_curr * m.R[2] + m.c_curr * m.R[3]
        m.R[2] = m.c_curr * m.R[2] + m.s_curr * m.R[3]
        m.R[3] = tmp
    end

    # Next rotation
    c, s, m.R[3] = givensAlgorithm(m.R[3], m.R[4])

    # Apply as well to the right-hand side
    m.rhs[2] = -s * m.rhs[1]
    m.rhs[1] = c * m.rhs[1]

    # Update W = V * inv(R). Two axpy's can maybe be one MV.
    copy!(m.w_next, m.v_curr)
    iteration > 1 && axpy!(-m.R[2], m.w_curr, m.w_next)
    iteration > 2 && axpy!(-m.R[1], m.w_prev, m.w_next)
    scale!(m.w_next, inv(m.R[3]))

    # Update solution x
    axpy!(m.rhs[1], m.w_next, m.x)

    # Move on: next -> curr, curr -> prev
    m.v_prev, m.v_curr, m.v_next = m.v_curr, m.v_next, m.v_prev
    m.w_prev, m.w_curr, m.w_next = m.w_curr, m.w_next, m.w_prev
    m.c_prev, m.s_prev, m.c_curr, m.s_curr = m.c_curr, m.s_curr, c, s
    m.rhs[1] = m.rhs[2]

    # Due to symmetry of the tri-diagonal matrix
    m.R[2] = m.prev_norm

    # The approximate residual is cheaply available
    m.resnorm = abs(m.rhs[2])

    m.resnorm, iteration + 1
end

function minres!(x, A, b; kwargs...)
    iterable = minres_iterable!(x, A, b; kwargs...)

    for resnorm = iterable end

    iterable.x, iterable.resnorm
end

minres(A, b; kwargs...) = minres!(zeros(b), A, b; initially_zero = true, kwargs...)