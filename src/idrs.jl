import Base.LinAlg.BLAS: axpy!
export idrs, idrs!

####

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

"""
IDRsSolver.jl
============================

The Induced Dimension Reduction method is a family of simple and fast Krylov
subspace algorithms for solving large nonsymmetric linear systems. The idea
behind the IDR(s) variant is to generate residuals that are in the nested
subspaces of shrinking dimension s.



Syntax
------

        X, h::ConvergenceHistory = idrs(A, b; s = 8, tol=sqrt(eps(typeof(real(b[1])))), maxiter = length(b)^2))

        The operator A must support multiplication with b,
        the right hand side b must support vecnorm, vecdot, copy, rand! and axpy!.
        .

Arguments
---------

       s -- dimension reduction number. Normally, a higher s gives faster convergence,
            but also  makes the method more expensive.
       tol -- tolerance of the method.
       maxiter -- maximum number of iterations


Output
------

        X -- Approximated solution by IDR(s)
        h -- Convergence history


        The [`ConvergenceHistory`](https://github.com/JuliaLang/IterativeSolvers.jl/issues/6) type provides information about the iteration history.
            - `isconverged::Bool`, a flag for whether or not the algorithm is converged.
            - `threshold`, the convergence threshold
            - `residuals::Vector`, the value of the convergence criteria at each iteration




References
----------

    [1] IDR(s): a family of simple and fast algorithms for solving large
        nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen
        SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, 2008
    [2] Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits
        Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld
        ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011
    [3] This file is a translation of the following MATLAB implementation:
        http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m
    [4] IDR(s)' webpage http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html
"""
idrs(A, b; s = 8, tol=sqrt(eps(typeof(real(b[1])))), maxiter = length(b)^2) =
    idrs_core!(zerox(A,b), linsys_op, (A,), b, s, tol, maxiter)

idrs!(x, A, b; s = 8, tol=sqrt(eps(typeof(real(b[1])))), maxiter=length(x)^2) =
    idrs_core!(x, linsys_op, (A,), b, s, tol, maxiter)

function idrs_core!{T}(X::T, op, args, C::T, s::Int64, tol::Float64, maxiter::Int64; smoothing::Bool=false)
    R = C - op(X, args...)::T
    normR = vecnorm(R)
    res = typeof(tol)[normR]
	iter = 0

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
    iter = 0
    while normR > tol && iter < maxiter
        for i in 1:s,
            f[i] = vecdot(P[i], R)
        end
        for k in 1:s

            # Solve small system and make v orthogonal to P

            c = LowerTriangular(M[k:s,k:s])\f[k:s]
            copy!(V, G[k])
            scale!(c[1], V)

            copy!(Q, U[k])
            scale!(c[1], Q)
            for i = k+1:s
                axpy!(c[i-k+1], G[i], V)
                axpy!(c[i-k+1], U[i], Q)
            end

            # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j

            #V = R - V
            scale!(-1., V)
            axpy!(1., R, V)

            copy!(U[k], Q)
            axpy!(om, V, U[k])
            G[k] = op(U[k], args...)

            # Bi-orthogonalise the new basis vectors

            for i in 1:k-1
                alpha = vecdot(P[i],G[k])/M[i,i]
                axpy!(-alpha, G[i], G[k])
                axpy!(-alpha, U[i], U[k])
            end

            # New column of M = P'*G  (first k-1 entries are zero)

            for i in k:s
                M[i,k] = vecdot(P[i],G[k])
            end

            #  Make r orthogonal to q_i, i = 1..k

            beta = f[k]/M[k,k]
            axpy!(-beta, G[k], R)
            axpy!(beta, U[k], X)

            normR = vecnorm(R)
            if smoothing
                copy!(T_s, R_s); axpy!(-1., R, T_s) # T_s = R_s - R
                gamma = vecdot(R_s, T_s)/vecnorm(T_s)^2
                axpy!(-gamma, T_s, R_s) # R_s = R_s - gamma*T_s
                copy!(T_s, X_s); axpy!(-1., X, T_s); axpy!(-gamma, T_s, X_s) # X_s = X_s - gamma*(X_s - X)
                normR = vecnorm(R_s)
            end
            res = [res; normR]
            iter += 1
            if normR < tol || iter > maxiter
                return X, ConvergenceHistory(normR < tol, tol, length(res), res)
            end
            if k < s
                f[k+1:s] = f[k+1:s] - beta*M[k+1:s,k]
            end

        end

        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r
        copy!(V, R)
        Q = op(V, args...)::T
        om = omega(Q, R)
        axpy!(-om, Q, R)
        axpy!(om, V, X)

        normR = vecnorm(R)
        if smoothing
            copy!(T_s, R_s); axpy!(-1., R, T_s) # T_s = R_s - R
            gamma = vecdot(R_s, T_s)/vecnorm(T_s)^2
            axpy!(-gamma, T_s, R_s) # R_s = R_s - gamma*T_s
            copy!(T_s, X_s); axpy!(-1., X, T_s); axpy!(-gamma, T_s, X_s) # X_s = X_s - gamma*(X_s - X)
            normR = vecnorm(R_s)
        end
        iter += 1
        res = [res; normR]
    end
    if smoothing
        X = copy(X_s)
    end
    return X, ConvergenceHistory(res[end]<tol, tol, length(res), res)
end

