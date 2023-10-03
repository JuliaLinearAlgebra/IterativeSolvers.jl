export idrs, idrs!
import Base: iterate
using Printf
using Random

"""
    idrs(A, b; s = 8, kwargs...) -> x, [history]

Same as [`idrs!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
idrs(A, b; kwargs...) = idrs!(zerox(A,b), A, b; kwargs...)

"""
    idrs!(x, A, b; s = 8, kwargs...) -> x, [history]

Solve the problem ``Ax = b`` approximately with IDR(s), where `s` is the dimension of the
shadow space.

# Arguments

- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `s::Integer = 8`: dimension of the shadow space;
- `Pl::precT`: left preconditioner,
- `abstol::Real = zero(real(eltype(b)))`,
  `reltol::Real = sqrt(eps(real(eltype(b))))`: absolute and relative
  tolerance for the stopping condition
  `|r_k| ≤ max(reltol * |r_0|, abstol)`, where `r_k = A * x_k - b`
  is the residual in the `k`th iteration;
- `maxiter::Int = size(A, 2)`: maximum number of iterations;
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iterations.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;
- `history`: convergence history.
"""
function idrs!(x, A, b;
               s = 8,
               Pl = Identity(),
               abstol::Real = zero(real(eltype(b))),
               reltol::Real = sqrt(eps(real(eltype(b)))),
               maxiter = size(A, 2),
               log::Bool = false,
               kwargs...)
    history = ConvergenceHistory(partial=!log)
    history[:abstol] = abstol
    history[:reltol] = reltol
    log && reserve!(history, :resnorm, maxiter)
    idrs_method!(history, x, A, b, s, Pl, abstol, reltol, maxiter; kwargs...)
    log && shrink!(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

@inline function omega(t, s)
    angle = sqrt(2.)/2
    ns = norm(s)
    nt = norm(t)
    ts = dot(t,s)
    rho = abs(ts/(nt*ns))
    omega = ts/(nt*nt)
    if rho < angle
        omega = omega*convert(typeof(omega),angle)/rho
    end
    omega
end

mutable struct IDRSIterable{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,
                            T15,T16,T17,T18,T19,T20,T21,T22,T23,T24,T25,T26}
    X::T1
    A::T2
    s::T3
    Pl::T4
    abstol::T5
    reltol::T6
    maxiter::T7
    smoothing::T8
    verbose::T9
    R::T10
    X_s::T11
    R_s::T12
    T_s::T13
    normR::T14
    tol::T15
    Z::T16
    P::T17
    U::T18
    G::T19
    Q::T20
    V::T21
    M::T22
    f::T23
    c::T24
    omega::T25
    log::T26
end
function idrs_iterable!(log, X, A, C::T,
    s::Number, Pl::precT, abstol::Real, reltol::Real, maxiter::Number; smoothing::Bool=false, verbose::Bool=false
    ) where {T, precT}
    R = C - A*X
    normR = norm(R)
    tol = max(reltol * normR, abstol)

    if smoothing
        X_s = copy(X)
        R_s = copy(R)
        T_s = zero(R)
    else 
        X_s = nothing
        R_s = nothing
        T_s = nothing
    end


    Z = zero(C)

    P = T[rand!(copy(C)) for k in 1:s]
    U = T[copy(Z) for k in 1:s]
    G = T[copy(Z) for k in 1:s]
    Q = copy(Z)
    V = copy(Z)

    M = Matrix{eltype(C)}(I,s,s)
    f = zeros(eltype(C),s)
    c = zeros(eltype(C),s)

    omega::eltype(C) = 1
    IDRSIterable(X, A, s, Pl, abstol, reltol, maxiter, smoothing, verbose,
    R, X_s, R_s, T_s, normR, tol, Z, P, U, G, Q, V, M, f, c, omega, log)
end


function idrs_method!(log::ConvergenceHistory, X, A, C::T,
    s::Number, Pl::precT, abstol::Real, reltol::Real, maxiter::Number; smoothing::Bool=false, verbose::Bool=false
    ) where {T, precT}

    verbose && @printf("=== idrs ===\n%4s\t%4s\t%7s\n", "iter", "step", "resnorm")

    iterable = idrs_iterable!(log, X, A, C, s, Pl, abstol, reltol, maxiter; smoothing, verbose)
    
    normR = reduce((_,r) -> r, iterable; init=iterable.normR)

    verbose && @printf("\n")
    iterable.X
end

function iterate(it::IDRSIterable, (iter, step) = (1, 1))
    X, A, s, Pl, R, X_s, R_s, T_s, Z, P, U, G, Q, V, M, f, c = 
    it.X, it.A, it.s, it.Pl, it.R, it.X_s, it.R_s, it.T_s, it.Z, it.P, it.U, it.G, it.Q, it.V, it.M, it.f, it.c
    
    if it.normR < it.tol || iter > it.maxiter
        it.log !== nothing && setconv(it.log, 0 <= it.normR < it.tol)

        if it.smoothing
            copyto!(X, X_s)
        end
        return nothing
    end

    if step in 1:s
        if step == 1     
            for i in 1:s
                f[i] = dot(P[i], R)
            end
        end
        k = step 

        # Solve small system and make v orthogonal to P

        c = LowerTriangular(M[k:s,k:s])\f[k:s]
        V .= c[1] .* G[k]
        Q .= c[1] .* U[k]

        for i = k+1:s
            V .+= c[i-k+1] .* G[i]
            Q .+= c[i-k+1] .* U[i]
        end

        # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j
        V .= R .- V

        # Preconditioning
        ldiv!(Pl, V)

        U[k] .= Q .+ it.omega .* V
        mul!(G[k], A, U[k])

        # Bi-orthogonalise the new basis vectors

        for i in 1:k-1
            alpha = dot(P[i], G[k])/M[i,i]
            G[k] .-= alpha .* G[i]
            U[k] .-= alpha .* U[i]
        end

        # New column of M = P'*G  (first k-1 entries are zero)

        for i in k:s
            M[i,k] = dot(P[i], G[k])
        end

        #  Make r orthogonal to q_i, i = 1..k

        beta = f[k]/M[k,k]
        R .-= beta .* G[k]
        X .+= beta .* U[k]

        it.normR = norm(R)
        if it.smoothing
            T_s .= R_s .- R

            gamma = dot(R_s, T_s)/dot(T_s, T_s)

            R_s .-= gamma .* T_s
            X_s .-= gamma .* (X_s .- X)

            it.normR = norm(R_s)
        end
        if k < s
            f[k+1:s] .-= beta*M[k+1:s,k]
        end
        nextstep = step + 1
    elseif step == s + 1

        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r
        copyto!(V, R)

        # Preconditioning
        ldiv!(Pl, V)

        mul!(Q, A, V)
        it.omega = omega(Q, R)
        R .-= it.omega .* Q
        X .+= it.omega .* V

        it.normR = norm(R)
        if it.smoothing
            T_s .= R_s .- R

            gamma = dot(R_s, T_s)/dot(T_s, T_s)

            R_s .-= gamma .* T_s
            X_s .-= gamma .* (X_s .- X)

            it.normR = norm(R_s)
        end
        nextstep = 1
    end
    if it.log !== nothing
        nextiter!(it.log, mvps=1)
        push!(it.log, :resnorm, it.normR)
    end
    it.verbose && @printf("%3d\t%3d\t%1.2e\n", iter, step, it.normR)
    return it.normR, (iter + 1, nextstep)
end

