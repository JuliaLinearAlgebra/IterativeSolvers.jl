export improved_gmres, improved_gmres!

using Base.LinAlg.Givens

struct NoopPreconditioner end

Base.A_ldiv_B!(::NoopPreconditioner, x) = x
Base.A_ldiv_B!(y, ::NoopPreconditioner, x) = copy!(y, x)

improved_gmres(A, b; kwargs...) = improved_gmres!(zeros(b), A, b; kwargs...)

function improved_gmres!(x, A, b;
  Pl = NoopPreconditioner(),
  Pr = NoopPreconditioner(),
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
  reserve!(history, :resnorm, outer * restart)
  improved_gmres_method!(history, x, A, b; tol = tol, outer = outer, restart = restart, kwargs...)
  (plot || log) && shrink!(history)
  plot && showplot(history)
  log ? (x, history) : x
end

function improved_gmres_method!(log::ConvergenceHistory, x, A, b;
  Pl = NoopPreconditioner(),
  Pr = NoopPreconditioner(),
  tol = sqrt(eps(real(eltype(b)))),
  restart::Int = min(20, length(b)),
  outer::Int = 1,
  maxiter::Int = restart,
  verbose::Bool = false
)
    T = eltype(b)

    # Approximate solution
    arnoldi = ArnoldiDecomp(A, restart, T)
    residual = Residual(restart, T)

    reserved_vec = similar(b)

    # Initial residual
    r0 = b - A * x

    # Preconditioned initial residual
    A_ldiv_B!(Pl, r0)

    # Log the first mvp
    log.mvps += 1

    residual.current = norm(r0)

    # Stopping criterion is based on |r0| / |rk|
    reltol = residual.current * tol

    for restarts = 1 : outer
    
        if residual.current ≤ reltol
          setconv(log, true)
          break
        end
        
        # Set the first basis vector
        β::real(T) = init!(arnoldi, x, b, Pl, reserved_vec)

        log.mvps += 1

        # And initialize the residual
        init_residual!(residual, β)

        # Inner iterations k = 1, ..., restart
        k = 1

        while k ≤ restart && residual.current > reltol

            # Arnoldi step: expand and orthogonalize
            expand!(arnoldi, Pl, Pr, k)
            orthogonalize!(arnoldi, k)
            update_residual!(residual, arnoldi, k)
            
            # Logging
            nextiter!(log, mvps = 1)
            push!(log, :resnorm, residual.current)
            verbose && @printf("%3d\t%3d\t%1.2e\n", restarts, k, residual.current)

            k += 1
        end

        # Solve the projected problem Hy = β * e1 in the least-squares sense.
        extract!(x, arnoldi, β, Pl, k)
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
  history::Vector{resT} # Residual per iteration
  current::resT # Current relative residual
  accumulator::resT # Used to compute the residual on the go
  nullvec::Vector{numT} # Vector in the null space to compute residuals
  β::resT # the initial residual
end

Residual(order, T::Type) = Residual{T, real(T)}(
  real(T)[],
  one(real(T)),
  one(real(T)),
  ones(T, order + 1),
  one(real(T))
)

function update_residual!(r::Residual, arnoldi::ArnoldiDecomp, k::Int)
  # Cheaply computes the current residual
  r.nullvec[k + 1] = -conj(dot(r.nullvec[1 : k], @view(arnoldi.H[1 : k, k])) / arnoldi.H[k + 1, k])
  r.accumulator += abs2(r.nullvec[k + 1])
  r.current = r.β / √r.accumulator
  push!(r.history, r.current)
end

function init!{T}(arnoldi::ArnoldiDecomp{T}, x, b, Pl, reserved_vec)
  # Initialize the Krylov subspace with the initial residual vector
  
  # This basically does V[1] = Pl \ (b - A * x)
  copy!(arnoldi.V[1], b)
  A_mul_B!(reserved_vec, arnoldi.A, x)
  axpy!(-one(T), reserved_vec, arnoldi.V[1])
  A_ldiv_B!(Pl, arnoldi.V[1])
  β = norm(arnoldi.V[1])
  scale!(arnoldi.V[1], one(T) / β)
  β
end

@inline function init_residual!{numT, resT}(r::Residual{numT, resT}, β)
  r.accumulator = one(resT)
  r.β = β
end

function extract!{T}(x, arnoldi::ArnoldiDecomp{T}, β, Pl, k::Int)
  # Computes & updates the solution
  rhs = zeros(T, k)
  rhs[1] = β
  y = solve!(Hessenberg(@view(arnoldi.H[1 : k, 1 : k - 1])), rhs)

  # Update x ← x + V * y
  for i = 1 : k - 1
    axpy!(y[i], arnoldi.V[i], x)
  end

  # Update x ← Pl \ x
  A_ldiv_B!(Pl, x)
end

@inline function expand!(arnoldi::ArnoldiDecomp, Pl, Pr, k::Int)
  # Simply expands by Pl * A * Pr * v
  A_mul_B!(arnoldi.V[k + 1], arnoldi.A, arnoldi.V[k])
  # A_ldiv_B!(arnoldi.V[k + 1], Pr, arnoldi.V[k])
  # copy!(arnoldi.V[k + 1], arnoldi.A * arnoldi.V[k + 1])
  A_ldiv_B!(Pl, arnoldi.V[k + 1])
end

function orthogonalize!{T}(arnoldi::ArnoldiDecomp{T}, k::Int)
  # Orthogonalize using Gram-Schmidt
  for j = 1 : k
    arnoldi.H[j, k] = dot(arnoldi.V[j], arnoldi.V[k + 1])
    axpy!(-arnoldi.H[j, k], arnoldi.V[j], arnoldi.V[k + 1])
  end

  arnoldi.H[k + 1, k] = norm(arnoldi.V[k + 1])
  scale!(arnoldi.V[k + 1], one(T) / arnoldi.H[k + 1, k])
end

type Hessenberg{T<:AbstractMatrix}
    H::T
end

@inline Base.size(H::Hessenberg, args...) = size(H.H, args...)

function Base.A_mul_B!(G::Givens, H::Hessenberg)
    m, n = size(H)
    @inbounds for i = G.i1 : n
        a1, a2 = H.H[G.i1, i], H.H[G.i2, i]
        H.H[G.i1, i] =       G.c  * a1 + G.s * a2
        H.H[G.i2, i] = -conj(G.s) * a1 + G.c * a2
    end
    return H
end

function solve!(H::Hessenberg, rhs)
    width = size(H, 2)

    # Hessenberg -> UpperTriangular; also apply to r.h.s.
    @inbounds for i = 1 : width
        rotation = Givens(i, i + 1, H.H[i, i], H.H[i + 1, i])
        A_mul_B!(rotation, H)
        A_mul_B!(rotation, rhs)
    end

    # Solve the upper triangular problem.
    U = UpperTriangular(@view(H.H[1 : width, 1 : width]))
    U \ @view(rhs[1 : width])
end