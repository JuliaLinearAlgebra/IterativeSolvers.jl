export idrs, idrs!

using Random

"""
    idrs(A, b; s = 8) -> x, [history]

Same as [`idrs!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
idrs(A, b; kwargs...) = idrs!(zerox(A,b), A, b; kwargs...)

"""
    idrs!(x, A, b; s = 8) -> x, [history]

Solve the problem ``Ax = b`` approximately with IDR(s), where `s` is the dimension of the
shadow space.

# Arguments

- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `s::Integer = 8`: dimension of the shadow space;
- `tol`: relative tolerance;
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
    s = 8, tol=sqrt(eps(real(eltype(b)))), maxiter=size(A, 2),
    log::Bool=false, kwargs...
    )
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    idrs_method!(history, x, A, b, s, tol, maxiter; kwargs...)
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
    om = ts/(nt*nt)
    if rho < angle
        om = om*convert(typeof(om),angle)/rho
    end
    om
end

function idrs_method!(log::ConvergenceHistory, X, A, C::T,
    s::Number, tol::Number, maxiter::Number; smoothing::Bool=false, verbose::Bool=false
    ) where {T}

    verbose && @printf("=== idrs ===\n%4s\t%7s\n","iter","resnorm")
    R = C - A*X
    normR = norm(R)
	iter = 1

    if smoothing
        X_s = copy(X)
        R_s = copy(R)
        T_s = zero(R)
    end

    if normR <= tol # Initial guess is a good enough solution
        return X
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

    om::eltype(C) = 1
    while normR > tol && iter ≤ maxiter
        for i in 1:s,
            f[i] = dot(P[i], R)
        end
        for k in 1:s
            nextiter!(log,mvps=1)

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

            U[k] .= Q .+ om .* V
            mul!(G[k], A, U[k])

            # Bi-orthogonalise the new basis vectors

            for i in 1:k-1
                alpha = dot(P[i],G[k])/M[i,i]
                G[k] .-= alpha .* G[i]
                U[k] .-= alpha .* U[i]
            end

            # New column of M = P'*G  (first k-1 entries are zero)

            for i in k:s
                M[i,k] = dot(P[i],G[k])
            end

            #  Make r orthogonal to q_i, i = 1..k

            beta = f[k]/M[k,k]
            R .-= beta .* G[k]
            X .+= beta .* U[k]

            normR = norm(R)
            if smoothing
                T_s .= R_s .- R

                gamma = dot(R_s, T_s)/dot(T_s, T_s)

                R_s .-= gamma .* T_s
                X_s .-= gamma .* (X_s .- X)

                normR = norm(R_s)
            end
            push!(log, :resnorm, normR)
            verbose && @printf("%3d\t%1.2e\n",iter,normR)
            if normR < tol || iter == maxiter
                shrink!(log)
                setconv(log, 0<=normR<tol)
                return X
            end
            if k < s
                f[k+1:s] .-=  beta*M[k+1:s,k]
            end
            iter += 1
        end

        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r
        copyto!(V, R)
        mul!(Q, A, V)
        om = omega(Q, R)
        R .-= om .* Q
        X .+= om .* V

        normR = norm(R)
        if smoothing
            T_s .= R_s .- R

            gamma = dot(R_s, T_s)/dot(T_s, T_s)

            R_s .-= gamma .* T_s
            X_s .-= gamma .* (X_s .- X)

            normR = norm(R_s)
        end
        iter += 1
        nextiter!(log, mvps=1)
        push!(log, :resnorm, normR)
    end
    if smoothing
        copyto!(X, X_s)
    end
    verbose && @printf("\n")
    setconv(log, 0<=normR<tol)
    X
end
