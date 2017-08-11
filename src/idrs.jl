export idrs, idrs!

"""
    idrs(A, b; s = 8) -> x, [history]

Same as [`idrs!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
idrs(A, b; kwargs...) = idrs!(zerox(A,b), A, b; kwargs...)

"""
    idrs(A, b; s = 8) -> x, [history]

Solve the problem ``Ax = b`` approximate with IDR(s), where `s` is the dimension of the
shadow space.

# Arguments

- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.

## Keywords

- `s::Integer = 8`: dimension of the shadow space;
- `tol`: relative tolerance;
- `maxiter::Int`: maximum number of iterations;
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
    s = 8, tol=sqrt(eps(real(eltype(b)))), maxiter=length(x)^2,
    log::Bool=false, kwargs...
    )
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    idrs_method!(history, x, linsys_op, (A,), b, s, tol, maxiter; kwargs...)
    log && shrink!(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

@inline function omega(t, s)
    angle = sqrt(2.)/2
    ns = vecnorm(s)
    nt = vecnorm(t)
    ts = vecdot(t,s)
    rho = abs(ts/(nt*ns))
    om = ts/(nt*nt)
    if rho < angle
        om = om*convert(typeof(om),angle)/rho
    end
    om
end

@inline linsys_op(x, A) = A*x

function idrs_method!(log::ConvergenceHistory, X, op, args, C::T,
    s::Number, tol::Number, maxiter::Number; smoothing::Bool=false, verbose::Bool=false
    ) where {T}

    verbose && @printf("=== idrs ===\n%4s\t%7s\n","iter","resnorm")
    R = C - op(X, args...)::T
    normR = vecnorm(R)
	iter = 1

    if smoothing
        X_s = copy(X)
        R_s = copy(R)
        T_s = zeros(R)
    end

    if normR <= tol           # Initial guess is a good enough solution
        return X, ConvergenceHistory(0<= res[end] < tol, tol, length(res), res)
    end

    Z = zero(C)

    P = T[rand!(copy(C)) for k in 1:s]
    U = T[copy(Z) for k in 1:s]
    G = T[copy(Z) for k in 1:s]
    Q = copy(Z)
    V = copy(Z)

    M = eye(eltype(C),s,s)
    f = zeros(eltype(C),s)
    c = zeros(eltype(C),s)

    om::eltype(C) = 1
    while normR > tol && iter ≤ maxiter
        for i in 1:s,
            f[i] = vecdot(P[i], R)
        end
        for k in 1:s
            nextiter!(log,mvps=1)

            # Solve small system and make v orthogonal to P

            c = LowerTriangular(M[k:s,k:s])\f[k:s]
            @blas! V = G[k]
            @blas! V *= c[1]

            @blas! Q = U[k]
            @blas! Q *= c[1]
            for i = k+1:s
                @blas! V += c[i-k+1]*G[i]
                @blas! Q += c[i-k+1]*U[i]
            end

            # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j

            #V = R - V
            @blas! V *= -1.
            @blas! V += R

            @blas! U[k] = Q
            @blas! U[k] += om*V
            G[k] = op(U[k], args...)

            # Bi-orthogonalise the new basis vectors

            for i in 1:k-1
                alpha = vecdot(P[i],G[k])/M[i,i]
                @blas! G[k] += -alpha*G[i]
                @blas! U[k] += -alpha*U[i]
            end

            # New column of M = P'*G  (first k-1 entries are zero)

            for i in k:s
                M[i,k] = vecdot(P[i],G[k])
            end

            #  Make r orthogonal to q_i, i = 1..k

            beta = f[k]/M[k,k]
            @blas! R += -beta*G[k]
            @blas! X += beta*U[k]

            normR = vecnorm(R)
            if smoothing
                # T_s = R_s - R
                @blas! T_s = R_s
                @blas! T_s += (-1.)*R

                gamma = vecdot(R_s, T_s)/vecdot(T_s, T_s)

                @blas! R_s += -gamma*T_s
                # X_s = X_s - gamma*(X_s - X)
                @blas! T_s = X_s
                @blas! T_s += (-1.)*X
                @blas! X_s += -gamma*T_s

                normR = vecnorm(R_s)
            end
            push!(log, :resnorm, normR)
            verbose && @printf("%3d\t%1.2e\n",iter,normR)
            if normR < tol || iter == maxiter
                shrink!(log)
                setconv(log, 0<=normR<tol)
                return X
            end
            if k < s
                f[k+1:s] = f[k+1:s] - beta*M[k+1:s,k]
            end
            iter += 1
        end

        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r
        @blas! V = R
        Q = op(V, args...)::T
        om = omega(Q, R)
        @blas! R += -om*Q
        @blas! X += om*V

        normR = vecnorm(R)
        if smoothing
            # T_s = R_s - R
            @blas! T_s = R_s
            @blas! T_s += (-1.)*R

            gamma = vecdot(R_s, T_s)/vecdot(T_s, T_s)

            @blas! R_s += -gamma*T_s
            # X_s = X_s - gamma*(X_s - X)
            @blas! T_s = X_s
            @blas! T_s += (-1.)*X
            @blas! X_s += -gamma*T_s

            normR = vecnorm(R_s)
        end
        iter += 1
        nextiter!(log,mvps=1)
        push!(log, :resnorm, normR)
    end
    if smoothing
        @blas! X = X_s
    end
    verbose && @printf("\n")
    setconv(log, 0<=normR<tol)
    X
end
