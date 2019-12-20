import Base: iterate

#Simple methods
export powm, powm!, invpowm, invpowm!

mutable struct PowerMethodIterable{matT, vecT <: AbstractVector, numT <: Number, eigvalT <: Number}
    A::matT
    x::vecT
    tol::numT
    maxiter::Int
    θ::eigvalT
    r::vecT
    Ax::vecT
    residual::numT
end


##
## Iterators
##

@inline converged(p::PowerMethodIterable) = p.residual ≤ p.tol

@inline start(p::PowerMethodIterable) = 0

@inline done(p::PowerMethodIterable, iteration::Int) = iteration > p.maxiter || converged(p)

function iterate(p::PowerMethodIterable, iteration::Int=start(p))
    if done(p, iteration) return nothing end

    mul!(p.Ax, p.A, p.x)

    # Rayleigh quotient θ = x'Ax
    p.θ = dot(p.x, p.Ax)

    # (Previous) residual vector r = Ax - λx
    copyto!(p.r, p.Ax)
    axpy!(-p.θ, p.x, p.r)

    # Normed residual
    p.residual = norm(p.r)

    # Normalize the next approximation
    copyto!(p.x, p.Ax)
    rmul!(p.x, one(eltype(p.x)) / norm(p.x))

    p.residual, iteration + 1
end

# Transforms the eigenvalue back whether shifted or inversed
@inline transform_eigenvalue(θ, inverse::Bool, σ) = σ + (inverse ? inv(θ) : θ)

function powm_iterable!(A, x; tol = eps(real(eltype(A))) * size(A, 2) ^ 3, maxiter = size(A, 1))
    T = eltype(x)
    PowerMethodIterable(A, x, tol, maxiter, zero(T), similar(x), similar(x), floatmax(real(T)))
end

"""
    powm(B; kwargs...) -> λ, x, [history]

See [`powm!`](@ref). Calls `powm!(B, x0; kwargs...)` with
`x0` initialized as a random, complex unit vector.
"""
function powm(B; kwargs...)
    x0 = rand(Complex{real(eltype(B))}, size(B, 1))
    rmul!(x0, one(eltype(B)) / norm(x0))
    powm!(B, x0; kwargs...)
end

"""
    powm!(B, x; shift = zero(eltype(B)), inverse::Bool = false, kwargs...) -> λ, x, [history]

By default finds the approximate eigenpair `(λ, x)` of `B` where `|λ|` is largest.

# Arguments
- `B`: linear map, see the note below.
- `x`: normalized initial guess. Don't forget to use complex arithmetic when necessary.

## Keywords
- `tol::Real = eps(real(eltype(B))) * size(B, 2) ^ 3`: stopping tolerance for the residual norm;
- `maxiter::Integer = size(B,2)`: maximum number of iterations;
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iterations.

!!! note "Shift-and-invert"
    When applying shift-and-invert to ``Ax = λx`` with `invert = true` and `shift = ...`, note
    that the role of `B * b` becomes computing `inv(A - shift I) * b`. So rather than
    passing the linear map ``A`` itself, pass a linear map `B` that has the action of
    shift-and-invert. The eigenvalue is transformed back to an eigenvalue of the actual
    matrix ``A``.

# Return values

**if `log` is `false`**
- `λ::Number` approximate eigenvalue computed as the Rayleigh quotient;
- `x::Vector` approximate eigenvector.

**if `log` is `true`**
- `λ::Number`: approximate eigenvalue computed as the Rayleigh quotient;
- `x::Vector`: approximate eigenvector;
- `history`: convergence history.

**ConvergenceHistory keys**

- `:tol` => `::Real`: stopping tolerance;
- `:resnom` => `::Vector`: residual norm at each iteration.

# Examples

```julia
using LinearMaps
σ = 1.0 + 1.3im
A = rand(ComplexF64, 50, 50)
F = lu(A - σ * I)
Fmap = LinearMap{ComplexF64}((y, x) -> ldiv!(y, F, x), 50, ismutating = true)
λ, x = powm(Fmap, inverse = true, shift = σ, tol = 1e-4, maxiter = 200)
```
"""
function powm!(B, x;
    tol = eps(real(eltype(B))) * size(B, 2) ^ 3,
    maxiter = size(B, 1),
    shift = zero(eltype(B)),
    inverse::Bool = false,
    log::Bool = false,
    verbose::Bool = false
)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    reserve!(history, :resnorm, maxiter)
    verbose && @printf("=== powm ===\n%4s\t%7s\n", "iter", "resnorm")

    iterable = powm_iterable!(B, x, tol = tol, maxiter = maxiter)

    for (iteration, residual) = enumerate(iterable)
        nextiter!(history, mvps = 1)
        verbose && @printf("%3d\t%1.2e\n", iteration, residual)
    end

    setconv(history, converged(iterable))

    verbose && println()

    log && shrink!(history)

    λ = transform_eigenvalue(iterable.θ, inverse, shift)
    x = iterable.x

    log ? (λ, x, history) : (λ, x)
end

"""
    invpowm(B; shift = σ, kwargs...) -> λ, x, [history]

Find the approximate eigenpair `(λ, x)` of ``A`` near `shift`, where `B`
is a linear map that has the effect `B * v = inv(A - σI) * v`.

The method calls `powm!(B, x0; inverse = true, shift = σ)` with
`x0` a random, complex unit vector. See [`powm!`](@ref)

# Examples

```julia
using LinearMaps
σ = 1.0 + 1.3im
A = rand(ComplexF64, 50, 50)
F = lu(A - σ * I)
Fmap = LinearMap{ComplexF64}((y, x) -> ldiv!(y, F, x), 50, ismutating = true)
λ, x = invpowm(Fmap, shift = σ, tol = 1e-4, maxiter = 200)
```
"""
function invpowm(B; kwargs...)
    x0 = rand(Complex{real(eltype(B))}, size(B, 1))
    rmul!(x0, one(eltype(B)) / norm(x0))
    invpowm!(B, x0; kwargs...)
end

"""
    invpowm!(B, x0; shift = σ, kwargs...) -> λ, x, [history]

Find the approximate eigenpair `(λ, x)` of ``A`` near `shift`, where `B`
is a linear map that has the effect `B * v = inv(A - σI) * v`.

The method calls `powm!(B, x0; inverse = true, shift = σ)`. See [`powm!`](@ref).
"""
invpowm!(B, x0; kwargs...) = powm!(B, x0; inverse = true, kwargs...)
