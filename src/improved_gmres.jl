export improved_gmres, improved_gmres!

improved_gmres(A, b; kwargs...) = improved_gmres!(zeros(b), A, b; kwargs...)

function improved_gmres!(x, A, b;
  Pl = Identity(),
  Pr = Identity(),
  outer::Int = 1,
  tol = sqrt(eps(real(eltype(b)))),
  restart::Int = min(20, length(b)),
  plot::Bool = false,
  log::Bool = false,
  kwargs...
)
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial = !log, restart = restart)
    history[:tol] = tol
    log && reserve!(history, :resnorm, outer * restart)
    improved_gmres_method!(history, x, A, b; Pl = Pl, Pr = Pr, tol = tol, outer = outer, restart = restart, log = log, kwargs...)
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

    for restarts = 1 : outer
    
        # Converged?
        if residual.current ≤ reltol
            setconv(history, true)
            break
        end
        
        # Set the first basis vector
        if restarts > 1
            β = init!(arnoldi, x, b, Pl, reserved_vec)
            
            if log
                history.mvps += 1
            end

            # And initialize the residual
            init_residual!(residual, β)
        end

        # Inner iterations k = 1, ..., restart
        k = 1

        while k ≤ restart && residual.current > reltol

            # Arnoldi step: expand and orthogonalize
            expand!(arnoldi, Pl, Pr, k)
            orthogonalize!(arnoldi, k)
            update_residual!(residual, arnoldi, k)
            
            # Logging
            if log
                nextiter!(history, mvps = 1)
                push!(history, :resnorm, residual.current)
            end
            
            verbose && @printf("%3d\t%3d\t%1.2e\n", restarts, k, residual.current)

            k += 1
        end

        # Solve the projected problem Hy = β * e1 in the least-squares sense,
        y = solve_least_squares!(arnoldi, β, k)
        update_solution!(x, y, arnoldi, Pr, k)
    end

    verbose && @printf("\n")
    x
end

mutable struct ArnoldiDecomp{T}
    A
    V::Vector{Vector{T}} # Orthonormal basis vectors
    H::Matrix{T}         # Hessenberg matrix
end

ArnoldiDecomp(A, order::Int, T::Type) = ArnoldiDecomp{T}(
    A,
    [zeros(T, size(A, 1)) for i = 1 : order + 1],
    zeros(T, order + 1, order)
)

mutable struct Residual{numT, resT}
    current::resT # Current relative residual
    accumulator::resT # Used to compute the residual on the go
    nullvec::Vector{numT} # Vector in the null space to compute residuals
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
    r.nullvec[k + 1] = -conj(dot(UnsafeVectorView(r.nullvec, 1 : k), UnsafeVectorView(arnoldi.H, 1 : k, k)) / arnoldi.H[k + 1, k])
    r.accumulator += abs2(r.nullvec[k + 1])
    r.current = r.β / √r.accumulator
end

function init!(arnoldi::ArnoldiDecomp{T}, x, b, Pl, reserved_vec) where {T}
    # Initialize the Krylov subspace with the initial residual vector
    # This basically does V[1] = Pl \ (b - A * x) and then normalize
    copy!(arnoldi.V[1], b)
    A_mul_B!(reserved_vec, arnoldi.A, x)
    @blas! arnoldi.V[1] -= one(T) * reserved_vec
    A_ldiv_B!(Pl, arnoldi.V[1])

    # Normalize
    β = norm(arnoldi.V[1])
    @blas! arnoldi.V[1] *= one(T) / β
    β
end

@inline function init_residual!(r::Residual{numT, resT}, β) where {numT,resT}
    r.accumulator = one(resT)
    r.β = β
end

function solve_least_squares!(arnoldi::ArnoldiDecomp{T}, β, k::Int) where {T}
    # Computes & updates the solution
    rhs = zeros(T, k)
    rhs[1] = β

    H = MyHessenberg(view(arnoldi.H, 1 : k, 1 : k - 1))
    y = solve!(H, rhs)
end

function update_solution!(x, y, arnoldi::ArnoldiDecomp{T}, Pr::Identity, k::Int) where {T}
    # Update x ← x + V * y
    for i = 1 : k - 1
        @blas! x += y[i] * arnoldi.V[i]
    end
end

function update_solution!(x, y, arnoldi::ArnoldiDecomp{T}, Pr, k::Int) where {T}
    # Allocates a temporary while computing x ← x + Pr \ (V * y)
    tmp = zeros(x)
    for i = 1 : k - 1
        @blas! tmp += y[i] * arnoldi.V[i]
    end
    @blas! x += one(T) * (Pr \ tmp)
end

function expand!(arnoldi::ArnoldiDecomp, Pl::Identity, Pr::Identity, k::Int)
    # Simply expands by A * v
    A_mul_B!(arnoldi.V[k + 1], arnoldi.A, arnoldi.V[k])
end

function expand!(arnoldi::ArnoldiDecomp, Pl, Pr::Identity, k::Int)
    # Expands by Pl \ (A * v)
    A_mul_B!(arnoldi.V[k + 1], arnoldi.A, arnoldi.V[k])
    A_ldiv_B!(Pl, arnoldi.V[k + 1])
end

function expand!(arnoldi::ArnoldiDecomp, Pl, Pr, k::Int)
    # Expands by Pl \ (A * (Pr \ v))
    A_ldiv_B!(arnoldi.V[k + 1], Pr, arnoldi.V[k])
    copy!(arnoldi.V[k + 1], arnoldi.A * arnoldi.V[k + 1])
    A_ldiv_B!(Pl, arnoldi.V[k + 1])
end

function orthogonalize!(arnoldi::ArnoldiDecomp{T}, k::Int) where {T}
    # Orthogonalize using modified Gram-Schmidt
    for j = 1 : k
        arnoldi.H[j, k] = dot(arnoldi.V[j], arnoldi.V[k + 1])
        @blas! arnoldi.V[k + 1] -= arnoldi.H[j, k] * arnoldi.V[j]
    end

    arnoldi.H[k + 1, k] = norm(arnoldi.V[k + 1])
    @blas! arnoldi.V[k + 1] *= one(T) / arnoldi.H[k + 1, k]
end
