export improved_gmres, improved_gmres!

improved_gmres(A, b; kwargs...) = improved_gmres!(zeros(b), A, b; kwargs...)

function improved_gmres!(x, A, b;
  Pl = Identity(),
  Pr = Identity(),
  tol = sqrt(eps(real(eltype(b)))),
  restart::Int = min(20, length(b)),
  maxiter::Int = restart,
  plot::Bool = false,
  log::Bool = false,
  kwargs...
)
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial = !log, restart = restart)
    history[:tol] = tol
    log && reserve!(history, :resnorm, maxiter)
    improved_gmres_method!(history, x, A, b; Pl = Pl, Pr = Pr, tol = tol, maxiter = maxiter, restart = restart, log = log, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

function improved_gmres_method!(history::ConvergenceHistory, x, A, b;
    Pl = Identity(),
    Pr = Identity(),
    tol = sqrt(eps(real(eltype(b)))),
    restart::Int = min(20, length(b)),
    outer::Int = 1,
    maxiter::Int = restart,
    verbose::Bool = false,
    log = false
)
    T = eltype(b)

    # Approximate solution
    arnoldi = ArnoldiDecomp(A, restart, T)
    residual = Residual(restart, T)

    # Workspace vector to reduce the # allocs.
    reserved_vec = similar(b)
    β = residual.current = init!(arnoldi, x, b, Pl, reserved_vec)
    init_residual!(residual, β)

    # Log the first mvp for computing the initial residual
    if log
        history.mvps += 1
    end

    # Stopping criterion is based on |r0| / |rk|
    reltol = residual.current * tol

    # Total iterations (not reset after restart)
    total_iter = 1

    while total_iter ≤ maxiter

        # We already have the initial residual
        if total_iter > 1

            # Set the first basis vector
            β = init!(arnoldi, x, b, Pl, reserved_vec)

            # And initialize the residual
            init_residual!(residual, β)
            
            if log
                history.mvps += 1
            end
        end

        # Inner iterations k = 1, ..., restart
        k = 1

        while residual.current > reltol && k ≤ restart && total_iter ≤ maxiter

            # Arnoldi step: expand
            expand!(arnoldi, Pl, Pr, k)

            # Orthogonalize V[:, k + 1] w.r.t. V[:, 1 : k]
            arnoldi.H[k + 1, k] = orthogonalize_and_normalize!(
                view(arnoldi.V, :, 1 : k),
                view(arnoldi.V, :, k + 1),
                view(arnoldi.H, 1 : k, k)
            )

            # Implicitly computes the residual
            update_residual!(residual, arnoldi, k)
            
            if log
                nextiter!(history, mvps = 1)
                push!(history, :resnorm, residual.current)
            end
            
            verbose && @printf("%3d\t%3d\t%1.2e\n", mod(total_iter, restart), k, residual.current)

            k += 1
            total_iter += 1
        end

        # Solve the projected problem Hy = β * e1 in the least-squares sense
        rhs = solve_least_squares!(arnoldi, β, k)

        # And improve the solution x ← x + Pr \ (V * y)
        update_solution!(x, view(rhs, 1 : k - 1), arnoldi, Pr, k)
    
        # Converged?
        if residual.current ≤ reltol
            setconv(history, true)
            break
        end
    end

    verbose && @printf("\n")
    x
end

type ArnoldiDecomp{T}
    A
    V::Matrix{T} # Orthonormal basis vectors
    H::Matrix{T} # Hessenberg matrix
end

ArnoldiDecomp(A, order::Int, T::Type) = ArnoldiDecomp{T}(
    A,
    zeros(T, size(A, 1), order + 1),
    zeros(T, order + 1, order)
)

type Residual{numT, resT}
    current::resT # Current relative residual
    accumulator::resT # Used to compute the residual on the go
    nullvec::Vector{numT} # Vector in the null space of H to compute residuals
    β::resT # the initial residual
end

Residual(order, T::Type) = Residual{T, real(T)}(
    one(real(T)),
    one(real(T)),
    ones(T, order + 1),
    one(real(T))
)

function update_residual!(r::Residual, arnoldi::ArnoldiDecomp, k::Int)
    # Cheaply computes the current residual
    r.nullvec[k + 1] = -conj(dot(view(r.nullvec, 1 : k), view(arnoldi.H, 1 : k, k)) / arnoldi.H[k + 1, k])
    r.accumulator += abs2(r.nullvec[k + 1])
    r.current = r.β / √r.accumulator
end

function init!{T}(arnoldi::ArnoldiDecomp{T}, x, b, Pl, reserved_vec)
    # Initialize the Krylov subspace with the initial residual vector
    # This basically does V[1] = Pl \ (b - A * x) and then normalize
    
    first_col = view(arnoldi.V, :, 1)

    copy!(first_col, b)
    A_mul_B!(reserved_vec, arnoldi.A, x)
    @blas! first_col -= one(T) * reserved_vec
    A_ldiv_B!(Pl, first_col)

    # Normalize
    β = norm(first_col)
    @blas! first_col *= one(T) / β
    β
end

@inline function init_residual!{numT,resT}(r::Residual{numT, resT}, β)
    r.accumulator = one(resT)
    r.β = β
end

function solve_least_squares!{T}(arnoldi::ArnoldiDecomp{T}, β, k::Int)
    # Compute the least-squares solution to Hy = β e1 via Given's rotations
    rhs = zeros(T, k)
    rhs[1] = β

    H = Hessenberg(view(arnoldi.H, 1 : k, 1 : k - 1))
    solve!(H, rhs)

    rhs
end

function update_solution!{T}(x, y, arnoldi::ArnoldiDecomp{T}, Pr::Identity, k::Int)
    # Update x ← x + V * y

    # TODO: find the SugarBLAS alternative
    BLAS.gemv!('N', one(T), view(arnoldi.V, :, 1 : k - 1), y, one(T), x)
end

function update_solution!{T}(x, y, arnoldi::ArnoldiDecomp{T}, Pr, k::Int)
    # Allocates a temporary while computing x ← x + Pr \ (V * y)
    tmp = view(arnoldi.V, :, 1 : k - 1) * y
    @blas! x += one(T) * (Pr \ tmp)
end

function expand!(arnoldi::ArnoldiDecomp, Pl::Identity, Pr::Identity, k::Int)
    # Simply expands by A * v without allocating
    A_mul_B!(view(arnoldi.V, :, k + 1), arnoldi.A, view(arnoldi.V, :, k))
end

function expand!(arnoldi::ArnoldiDecomp, Pl, Pr::Identity, k::Int)
    # Expands by Pl \ (A * v) without allocating
    A_mul_B!(view(arnoldi.V, :, k + 1), arnoldi.A, view(arnoldi.V, :, k))
    A_ldiv_B!(Pl, view(arnoldi.V, :, k + 1))
end

function expand!(arnoldi::ArnoldiDecomp, Pl, Pr, k::Int)
    # Expands by Pl \ (A * (Pr \ v)). Allocates one vector.
    A_ldiv_B!(view(arnoldi.V, :, k + 1), Pr, view(arnoldi.V, :, k))
    copy!(view(arnoldi.V, :, k + 1), arnoldi.A * view(arnoldi.V, :, k + 1))
    A_ldiv_B!(Pl, view(arnoldi.V, :, k + 1))
end
