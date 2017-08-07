import Base: start, next, done

#Simple methods
export powm, invpowm

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

function next(p::PowerMethodIterable, iteration::Int)

    A_mul_B!(p.Ax, p.A, p.x)

    # Rayleigh quotient θ = x'Ax
    p.θ = dot(p.x, p.Ax)

    # (Previous) residual vector r = Ax - λx
    copy!(p.r, p.Ax)
    BLAS.axpy!(-p.θ, p.x, p.r)

    # Normed residual
    p.residual = norm(p.r)

    # Normalize the next approximation
    copy!(p.x, p.Ax)
    scale!(p.x, one(eltype(p.x)) / norm(p.x))

    p.residual, iteration + 1
end

# Transforms the eigenvalue back whether shifted or inversed
@inline transform_eigenvalue(θ, inverse::Bool, σ) = σ + (inverse ? inv(θ) : θ)

function powm_iterable!(A, x; tol = eps(real(eltype(A))) * size(A, 2) ^ 3, maxiter = size(A, 1))
    T = eltype(x)
    PowerMethodIterable(A, x, tol, maxiter, zero(T), similar(x), similar(x), realmax(real(T)))
end

function powm_iterable(A; kwargs...)
    x0 = rand(Complex{real(eltype(A))}, size(A, 1))
    scale!(x0, one(eltype(A)) / norm(x0))
    powm_iterable!(A, x0; kwargs...)
end

####################
# API method calls #
####################

function powm(A; kwargs...)
    x0 = rand(Complex{real(eltype(A))}, size(A, 1))
    scale!(x0, one(eltype(A)) / norm(x0))
    powm!(A, x0; kwargs...)
end

function powm!(A, x;
    tol = eps(real(eltype(A))) * size(A, 2) ^ 3,
    maxiter = size(A, 1),
    shift = zero(eltype(A)),
    inverse::Bool = false,
    log::Bool = false,
    verbose::Bool = false
)
    history = ConvergenceHistory(partial = !log)
    history[:tol] = tol
    reserve!(history, :resnorm, maxiter)
    verbose && @printf("=== powm ===\n%4s\t%7s\n", "iter", "resnorm")

    iterable = powm_iterable!(A, x, tol = tol, maxiter = maxiter)

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

function invpowm(B; kwargs...)
    x0 = rand(Complex{real(eltype(B))}, size(B, 1))
    scale!(x0, one(eltype(B)) / norm(x0))
    invpowm!(B, x0; kwargs...)
end

invpowm!(B, x0; kwargs...) = powm!(B, x0; inverse = true, kwargs...)

#################
# Documentation #
#################

let
#Initialize parameters
doc1_call = """    powm(A)
"""
doc2_call = """    invpowm(B)
"""
doc1_msg = """Find the largest eigenvalue `λ` (in absolute value) of `A` and its associated eigenvector `x`
using the power method.
"""
doc2_msg = """For an eigenvalue problem Ax = λx, find the closest eigenvalue in the complex plane to `shift`
together with its associated eigenvector `x`. The first argument `B` should be a mapping that has
the effect of B * x = inv(A - shift * I) * x and should support the `A_mul_B!` operation.

# Examples

```julia
using LinearMaps
σ = 1.0 + 1.3im
A = rand(Complex128, 50, 50)
F = lufact(A - UniformScaling(σ))
Fmap = LinearMap{Complex128}((y, x) -> A_ldiv_B!(y, F, x), 50, ismutating = true)
λ, x = invpowm(Fmap, shift = σ, tol = 1e-4, maxiter = 200)
```
"""
doc1_karg = ""
doc2_karg = "`shift::Number=0`: shift to be applied to matrix A."

doc1_version = (powm, doc1_call, doc1_msg, doc1_karg)
doc2_version = (invpowm, doc2_call, doc2_msg, doc2_karg)

i=0
docstring = Vector(2)

#Build docs
for (func, call, msg, karg) in [doc1_version, doc2_version]
i+=1
docstring[i] = """
$call

$msg

If `log` is set to `true` is given, method will output a tuple `λ, x, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `λ, x`.

# Arguments

`A`: linear operator.

## Keywords

$karg

`x = random unit vector`: initial eigenvector guess.

`tol::Real = eps()*size(A,2)^3`: stopping tolerance.

`maxiter::Integer = size(A,2)`: maximum number of iterations.

`verbose::Bool = false`: verbose flag.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

# Output

**if `log` is `false`**

`λ::Number`: eigenvalue

`x::Vector`: eigenvector

**if `log` is `true`**

`eig::Real`: eigenvalue

`x::Vector`: eigenvector

`ch`: convergence history.

**ConvergenceHistory keys**

`:tol` => `::Real`: stopping tolerance.

`:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> powm
@doc docstring[2] -> invpowm
end
