#################################################################
# Randomized estimators of elementary linear algebraic quantities
#################################################################

export rcond, reigmax, reigmin, rnorm, rnorms


"""
    randnn(el, m)
    randnn(el, m, n)

Compute randomized gaussian matrix normalized by column.

# Arguments

`el::Type`: element type.

`m::Int`: number of rows.

`n::Int`: number of columns.

## Keywords

`normalize::Bool = true`: normalize output.

# Output

**Without `n`**

`Ω`: vector containing Gaussian random numbers of type `el`.

**With `n`**

`Ω`: matrix of dimensions `m` x `n` containing Gaussian random numbers of
type `el`.
"""
function randnn(el, m::Int, normalize::Bool=true)
    if el <: Real
        Ω = randn(m)
    elseif el <: Complex
        Ω = randn(m) + im*randn(m)
    else
        throw(ArgumentError("Unsupported element type: $el"))
    end
    normalize ? Ω/norm(Ω) : Ω
end
function randnn(el, m::Int, n::Int, normalize::Bool=true)
    if el <: Real
        Ω = randn(m, n)
    elseif el <: Complex
        Ω = randn(m, n) + im*randn(m, n)
    else
        throw(ArgumentError("Unsupported element type: $el"))
    end
    normalize || return Ω
    for i=1:n
        Ω[:, i] /= norm(view(Ω, :, i))
    end
    Ω
end

"""
    rnorm(A, mvps)

Compute a probabilistic upper bound on the norm of a matrix `A`.
`‖A‖ ≤ α √(2/π) maxᵢ ‖Aωᵢ‖` with probability `p=α^(-mvps)`.

# Arguments

`A`: matrix whose norm to estimate.

`mvps::Int`: number of matrix-vector products to compute.

## Keywords

`p::Real=0.05`: probability of upper bound failing.

# Output

Estimate of ‖A‖.

# See also

see [`rnorms`](@ref) for a different estimator that uses premultiplying by both
`A` and `A'`.

# References

\cite[Lemma 4.1]{Halko2011}
"""
function rnorm(A, r::Int, p::Real=0.05)
    @assert 0<p≤1
    α = p^(-1.0/r)
    m, n = size(A)
    Ω = randnn(eltype(A), n, r, false)
    AΩ = A*Ω
    mx = maximum([norm(view(AΩ, :, j)) for j=1:r])
    α * √(2/π) * mx
end

"""
    rnorms(A, iters=1)

Estimate matrix norm randomly using `A'A`.

Compute a probabilistic upper bound on the norm of a matrix `A`.

	ρ = √(‖(A'A)ʲω‖/‖(A'A)ʲ⁻¹ω‖)

which is an estimate of the spectral norm of `A` produced by `iters`
steps of the power method starting with normalized `ω`, is a lower
bound on the true norm by a factor

	ρ ≤ α ‖A‖

with probability greater than `1 - p`, where `p = 4\sqrt(n/(iters-1)) α^(-2iters)`.

# Arguments

`A`: matrix whose norm to estimate.

`iters::Int = 1`: mumber of power iterations to perform.

## Keywords

`p::Real = 0.05`: probability of upper bound failing.

`At = A'`: Transpose of `A`.

# Output

Estimate of ‖A‖.

# See also

see [`rnorm`](@ref) for a different estimator that does not require
premultiplying by `A'`

# References

Appendix of \cite{Liberty2007}.

```bibtex
@article{Liberty2007,
    authors = {Edo Liberty and Franco Woolfe and Per-Gunnar Martinsson
    and Vladimir Rokhlin and Mark Tygert},
    title = {Randomized algorithms for the low-rank approximation of matrices},
    journal = {Proceedings of the National Academy of Sciences},
    volume = {104},
    issue = {51},
    year = {2007},
    pages = {20167--20172},
    doi  = {10.1073/pnas.0709640104}
}
```
"""
function rnorms(A, j::Int=1, p::Real=0.05; At = A')
    @assert 0<p≤1
     m, n = size(A)
    α = ((j-1)/n*(p/4)^2)^(-1/(4j))

    Ωold = Ω = randnn(eltype(A), n)
    for i=1:j #Power iterations
        Ω, Ωold = At*(A*Ω), Ω
    end
    ρ = √(norm(Ω)/norm(Ωold))
    α*ρ
end

"""
    rcond(A, iters=1)

Estimate matrix condition number randomly.

# Arguments

`A`: matrix whose condition number to estimate. Must be square and
support premultiply (`A*⋅`) and solve (`A\⋅`).

`iters::Int = 1`: number of power iterations to run.

## Keywords

`p::Real = 0.05`: probability that estimate fails to hold as an upper bound.

# Output

Interval `(x, y)` which contains `κ(A)` with probability `1 - p`.

# Implementation note

\cite{Dixon1983} originally describes this as a computation that
can be done by computing the necessary number of power iterations given p
and the desired accuracy parameter `θ=y/x`. However, these bounds were only
derived under the assumptions of exact arithmetic. Empirically, `iters≥4` has
been seen to result in incorrect results in that the computed interval does
not contain the true condition number. This implemention therefore makes `iters`
an explicitly user-controllable parameter from which to infer the accuracy
parameter and hence the interval containing `κ(A)`.

# References

\cite[Theorem 2]{Dixon1983}

```bibtex
@article{Dixon1983,
    author = {Dixon, John D},
    doi = {10.1137/0720053},
    journal = {SIAM Journal on Numerical Analysis},
    number = {4},
    pages = {812--814},
    title = {Estimating Extremal Eigenvalues and Condition Numbers of
	Matrices},
    volume = {20},
    year = {1983}
}
```
"""
function rcond(A, k::Int=1, p::Real=0.05)
    @assert 0<p≤1
    m, n = size(A)
    @assert m==n
    θ = (8n/(π*p^2))^(1/k)
    x = randnn(eltype(A), n)
    for i=1:k
        x = A*x
    end
    y = randnn(eltype(A), n)
    for i=1:k
        y = A\y
    end
    φ = ((x⋅x)*(y⋅y))^(1/(2k))
    (φ, θ*φ)
end

"""
    reigmax(A, iters=1)

Estimate maximal eigenvalue randomly.

# Arguments

`A`: Matrix whose maximal eigenvalue to estimate.
Must be square and support premultiply (`A*⋅`).

`iters::Int=1`: Number of power iterations to run. (Recommended: `iters` ≤ 3)

## Keywords

`p::Real=0.05`: Probability that estimate fails to hold as an upper bound.

# Output

Interval `(x, y)` which contains the maximal eigenvalue of `A` with
probability `1 - p`.

# References

\cite[Corollary of Theorem 1]{Dixon1983}.
"""
function reigmax(A, k::Int=1, p::Real=0.05)
    @assert 0<p≤1
    m, n = size(A)
    @assert m==n
    θ = (2n/(π*p^2))^(1/k)
    y = x = randnn(eltype(A), n)
    for i=1:k
        x = A*x
    end
    φ = y⋅x
    (φ, θ*φ)
end

"""
    reigmin(A, iters=1)

Estimate minimal eigenvalue randomly.

# Arguments

`A`: Matrix whose maximal eigenvalue to estimate.
Must be square and support premultiply (`A*⋅`).

`iters::Int=1`: Number of power iterations to run. (Recommended: `iters` ≤ 3)

## Keywords

`p::Real=0.05`: Probability that estimate fails to hold as an upper bound.

# Output

Interval `(x, y)` which contains the maximal eigenvalue of `A` with
probability `1 - p`.

# References

\cite[Corollary of Theorem 1]{Dixon1983}.
"""
function reigmin(A, k::Int=1, p::Real=0.05)
    @assert 0<p≤1
    m, n = size(A)
    @assert m==n
    θ = (2n/(π*p^2))^(1/k)
    y = x = randnn(eltype(A), n)
    for i=1:k
        x = A\x
    end
    φ = y⋅x
    (inv(θ*φ), inv(φ))
end

"""
    srft

A subsampled random Fourier transform.

# Fields

`l<:Integer`: number of vectors to return.

# Implements

`Base`: `*`
"""
immutable srft{T<:Integer}
    l :: T
end

#(*)(A::LinearMap, B::srft) = error("method only defined to avoid ambiguity. If you need this method please open a pull request")

"""
    *

Apply a subsampled random Fourier transform to the columns of `A`.

# Arguments

`A`: matrix to transform.

`Ω::srft`: subsampled random Fourier transform.

# Output

`B`: A matrix of dimensions size(A,1) x Ω.l

# References

\[Equation 4.6]{Halko2011}
"""
#Define two methods here to avoid method ambiguity with f::Function*b::Any
*(A::Function, Ω::srft) = function *(A, Ω::srft)
    m, n = size(A)
    B = A*Diagonal(exp(2π*im*rand(n))/√Ω.l)
    B = vcat([fft(A[i,:]) for i=1:m]...) #Factor of √n cancels out
    B[:, randperm(n)[1:Ω.l]]
end
