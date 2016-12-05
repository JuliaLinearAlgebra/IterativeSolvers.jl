import Base.LinAlg.BLAS: axpy!
export idrs, idrs!

####################
# API method calls #
####################


idrs(A, b; kwargs...) = idrs!(zerox(A,b), A, b; kwargs...)

function idrs!(x, A, b;
    s = 8, tol=sqrt(eps(typeof(real(b[1])))), maxiter=length(x)^2,
    plot::Bool=false, log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    idrs_method!(history, x, linsys_op, (A,), b, s, tol, maxiter; kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
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

"""
The Induced Dimension Reduction method is a family of simple and fast Krylov
subspace algorithms for solving large nonsymmetric linear systems. The idea
behind the IDR(s) variant is to generate residuals that are in the nested
subspaces of shrinking dimensions.
"""
function idrs_method!{T}(log::ConvergenceHistory, X, op, args, C::T,
    s::Number, tol::Number, maxiter::Number; smoothing::Bool=false, verbose::Bool=false
    )

    verbose && @printf("=== idrs ===\n%4s\t%7s\n","iter","resnorm")
    R = C - op(X, args...)::T
    normR = vecnorm(R)
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
    while normR > tol && iter < maxiter
        for i in 1:s,
            f[i] = vecdot(P[i], R)
        end
        for k in 1:s
            nextiter!(log,mvps=1)

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
                # T_s = R_s - R
                copy!(T_s, R_s); axpy!(-1., R, T_s)

                gamma = vecdot(R_s, T_s)/vecdot(T_s, T_s)

                # R_s = R_s - gamma*T_s
                axpy!(-gamma, T_s, R_s)
                # X_s = X_s - gamma*(X_s - X)
                copy!(T_s, X_s); axpy!(-1., X, T_s); axpy!(-gamma, T_s, X_s)

                normR = vecnorm(R_s)
            end
            push!(log, :resnorm, normR)
            verbose && @printf("%3d\t%1.2e\n",iter,normR)
            iter += 1
            if normR < tol || iter > maxiter
                shrink!(log)
                setconv(log, 0<=normR<tol)
                return X
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
            # T_s = R_s - R
            copy!(T_s, R_s); axpy!(-1., R, T_s)

            gamma = vecdot(R_s, T_s)/vecdot(T_s, T_s)

            # R_s = R_s - gamma*T_s
            axpy!(-gamma, T_s, R_s)
            # X_s = X_s - gamma*(X_s - X)
            copy!(T_s, X_s); axpy!(-1., X, T_s); axpy!(-gamma, T_s, X_s)

            normR = vecnorm(R_s)
        end
        iter += 1
        nextiter!(log,mvps=1)
        push!(log, :resnorm, normR)
    end
    if smoothing
        copy!(X, X_s)
    end
    verbose && @printf("\n")
    setconv(log, 0<=normR<tol)
    X
end

#################
# Documentation #
#################

let
#Initialize parameters
doc_call = """    idrs(A, b)
"""
doc!_call = """    idrs!(x, A, b)
"""

doc_msg = "Solve A*x=b using the induced dimension reduction method."
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """* `x`: initial guess, overwrite final estimation."""

doc_version = (idrs, doc_call, doc_msg, doc_arg)
doc!_version = (idrs!, doc!_call, doc!_msg, doc!_arg)

i=0
docstring = Vector(2)

#Build docs
for (func, call, msg, arg) in [doc_version, doc!_version]
i+=1
docstring[i] =  """
$call
$msg
If `log` is set to `true` is given, method will output a tuple `x, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `x`.
The `plot` attribute can only be used when `log` is set version.
**Arguments**
$arg
* `A`: linear operator.
* `b`: right hand side.
*Keywords*
* `Pl = 1`: left preconditioner of the method.
* `Pr = 1`: left preconditioner of the method.
* `tol::Real = sqrt(eps())`: stopping tolerance.
* `restart::Integer = min(20,length(b))`: maximum number of iterations per restart.
* `maxiter::Integer = min(20,length(b))`: maximum number of iterations.
* `verbose::Bool = false`: print method information.
* `log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.
* `plot::Bool = false`: plot data. (Only when `log` is set)
**Output**
*`log` is `false`:*
* `x`: approximated solution.
*`log` is `true`:*
* `x`: approximated solution.
* `ch`: convergence history.
*ConvergenceHistory keys*
* `:tol` => `::Real`: stopping tolerance.
* `:resnom` => `::Vector`: residual norm at each iteration.
"""
end

@doc docstring[1] -> idrs
@doc docstring[2] -> idrs!
end
