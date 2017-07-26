#########################################################
## Randomized singular value and eigen- decompositions ##
########################################################

# This file provides a rudimentary implementation of the randomized singular
# value decomposition and spectral (eigen-) decomposition as described in
# \cite{Halko2011}.

import Base.LinAlg: Eigen, SVD

export rsvdfact, reig

"""
    rsvdfact(A, n, p=0)

Compute partial singular value decomposition of `A` using a randomized
algorithm.

# Arguments

`A`: input matrix.

`n::Int`: number of singular value/vector pairs to find.

`p::Int=0`: number of extra vectors to include in computation.

# Output

`::SVD`: singular value decomposition.

# Warning

This variant of the randomized singular value decomposition is the most
commonly found implementation but is not recommended for accurate
computations, as it often has trouble finding the `n` largest singular pairs,
but rather finds `n` large singular pairs which may not necessarily be the
largest.

# Implementation note

This function calls `rrange`, which uses naive randomized rangefinding to
compute a basis for a subspace of dimension `n` (Algorithm 4.1 of
\\cite{Halko2011}), followed by `svdfact_restricted()`, which computes the
exact SVD factorization on the restriction of `A` to this randomly selected
subspace (Algorithm 5.1 of \\cite{Halko2011}).

Alternatively, you can mix and match your own randomized algorithm using
any of the randomized range finding algorithms to find a suitable subspace
and feeding the result to one of the routines that computes the `SVD`
restricted to that subspace.

# References

```bibtex
@article{Halko2011,
    author = {Halko, N and Martinsson, P G and Tropp, J A},
    doi = {10.1137/090771806},
    journal = {SIAM Review},
    month = jan,
    number = {2},
    pages = {217--288},
    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},
    volume = {53},
    year = {2011}
}
```
"""
function rsvdfact(A, n::Int, p::Int=0)
    Q = rrange(A, n+p)
    svdfact_restricted(A, Q, n)
end

"""
    rsvdvals(A, n, p=0)

Compute partial singular value decomposition of `A` using a randomized
algorithm.

# Arguments

`A`: input matrix.

`n::Int`: number of singular value/vector pairs to find.

`p::Int=0`: number of extra vectors to include in computation.

# Output

`::Vector`: estimated singular values of `A`.

# Warning

This variant of the randomized singular value decomposition is the most
commonly found implementation but is not recommended for accurate
computations, as it often has trouble finding the `n` largest singular pairs,
but rather finds `n` large singular pairs which may not necessarily be the
largest.

# Implementation note

This function calls `rrange`, which uses naive randomized rangefinding to
compute a basis for a subspace of dimension `n` (Algorithm 4.1 of
\\cite{Halko2011}), followed by `svdfact_restricted()`, which computes the
exact SVD factorization on the restriction of `A` to this randomly selected
subspace (Algorithm 5.1 of \\cite{Halko2011}).

Alternatively, you can mix and match your own randomized algorithm using
any of the randomized range finding algorithms to find a suitable subspace
and feeding the result to one of the routines that computes the `SVD`
restricted to that subspace.

# References

```bibtex
@article{Halko2011,
    author = {Halko, N and Martinsson, P G and Tropp, J A},
    doi = {10.1137/090771806},
    journal = {SIAM Review},
    month = jan,
    number = {2},
    pages = {217--288},
    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},
    volume = {53},
    year = {2011}
}
```
"""
function rsvdvals(A, n::Int, p::Int=0)
    Q = rrange(A, n+p)
    svdvals_restricted(A, Q, n)
end

"""
    rrange(A, r=0)

Compute an orthonormal basis for a subspace of `A` of dimension `r` using
naive randomized rangefinding.

# Arguments

`A`: Input matrix. Must support `size(A)` and premultiply.

`r::Int = 0`: The number of basis vectors to compute.

# Output

`::Matrix`: matrix of dimension `size(A,1) x r` containing the basis
vectors of the computed subspace of `A`.

# Warning

The Reference explicitly discourages using this algorithm.

# Implementation note

Whereas \\cite{Halko2011} recommends classical Gram-Schmidt with double
reorthogonalization, we instead compute the basis with `qrfact()`, which
for dense `A` computes the QR factorization using Householder reflectors.
"""
function rrange(A, l::Int=0)
    m, n = size(A)
    if l > m
	    throw(ArgumentError("Cannot find $l linearly independent vectors of $m x $n matrix"))
    end
    Ω = randn(n, l)
    Y = A*Ω
    Q = full(qrfact!(Y)[:Q])
    @assert m==size(Q, 1)
    @assert l==size(Q, 2)
    return Q
end

"""
    rrange_adaptive(A, r, ϵ=eps(); maxiter=10)

Compute an orthogonal basis for a subspace of `A` of dimension `l` using adaptive
randomized rangefinding.

Similar to `rrange`, but determines the oversampling parameter adaptively given
a threshold `ϵ`.

# Arguments

`A`: input matrix. Must support `size(A)` and premultiply.

`r::Integer`: number of basis vectors to compute.

`ϵ::Real = eps()`: threshold to determine adaptive fitting.

## Keywords

`maxiter::Int = 10`: maximum number of iterations to run.

# Output

`::Matrix`: matrix of dimension `size(A,1)` x l containing the basis
vectors of the computed subspace of `A`.

# References

Algorithm 4.2 of \\cite{Halko2011}
"""
function rrange_adaptive(A, r::Integer, ϵ::Real=eps(); maxiter::Int=10)
    m, n = size(A)

    Ω = randn(n,r)
    Y = A*Ω
    Q = zeros(m,0)

    const tol=ϵ/(10*√(2/π))
    for j=1:maxiter
        Tol = maximum([norm(Y[:,i]) for i=j:(j+r-1)])
        Tol > tol || break

        y = view(Y,:,j)
        y = Y[:,j] = y - Q*(Q'y)
        q = y/norm(y)
        Q = [Q q]

        ω = randn(n)
        y = A*ω
        y = y - Q*(Q'y)
        Y = [Y y]
        Yb = view(Y, :, (j+1):(j+r-1))
        Yb = Y[:, (j+1):(j+r-1)] = Yb - q * q'Yb

        j==maxiter && warn("Maximum number of iterations reached with norm $Tol > $tol")
    end
    Q
end

"""
    rrange_si(A, l; At=A', q=0)

Compute an orthonormal basis for a subspace of `A` of dimension `l` using
randomized rangefinding by subspace iteration.

*Note:* Running with `q=0` is functionally equivalent to `rrange`.

# Arguments

`A`: input matrix. Must support `size(A)` and premultiply.

`l::Int`: number of basis vectors to compute.

## Keywords

`At = A'`: transpose of `A`.

`q::Int = 0`: number of subspace iterations.

`p` : oversampling parameter. The number of extra basis vectors to use
in the computation, which get discarded at the end.

# Output

`::Matrix`: A dense matrix of dimension `size(A,1) x l` containing the basis
vectors of the computed subspace of `A`.

# Implementation note

Whereas the Reference recommends classical Gram-Schmidt with double
reorthogonalization, we instead compute the basis with `qrfact()`, which
for dense A computes the QR factorization using Householder reflectors.

# References

Algorithm 4.4 of \\cite{Halko2011}
"""
function rrange_si(A, l::Int; At=A', q::Int=0)
    basis=x->full(qrfact(x)[:Q])
    n = size(A, 2)
    Ω = randn(n,l+p)
    Y = A*Ω
    Q = basis(Y)
    for iter=1:q #Do some power iterations
        Ỹ=At*Q
        Q̃=basis(Ỹ)
        Y=A*Q̃
        Q=basis(Y)
    end
    Q
end

"""
    rrange_f(A, l)

Compute an orthonormal basis for a subspace of `A` of dimension `l` using
naive randomized rangefinding using stochastic randomized Fourier transforms.

*Note:* similar to `rrange`, but does not use gaussian random matrices.

# Arguments

`A` : input matrix. Must support `size(A)` and premultiply.

`l::Int` : number of basis vectors to compute.

`p::Int = 0` : oversampling parameter. The number of extra basis vectors to use
in the computation, which get discarded at the end.

# Output

`::Matrix`: matrix of dimension `size(A,1)` x `l` containing the basis
vectors of the computed subspace of `A`.

# Implementation note

Whereas the Reference recommends classical Gram-Schmidt with double
reorthogonalization, we instead compute the basis with `qrfact()`, which
for dense `A` computes the QR factorization using Householder reflectors.

# References

Algorithm 4.5 of \\cite{Halko2011}
"""
function rrange_f(A, l::Int)
    n = size(A, 2)
    Ω = srft(l+p)
    Y = A*Ω
    Q = full(qrfact!(Y)[:Q])
end

"""
    svdfact_restricted(A, Q, n)

Compute the SVD factorization of `A` restricted to the subspace spanned by `Q`
using exact projection.

# Arguments

`A`: input matrix. Must support postmultiply.

`Q`: matrix containing basis vectors of the subspace whose restriction to is
desired.

# Output

`::SVD`: singular value decomposition.

# References

Algorithm 5.1 of \\cite{Halko2011}
"""
function svdfact_restricted(A, Q, n::Int)
    B=Q'A
    S=svdfact!(B)
    SVD((Q*S[:U])[:, 1:n], S[:S][1:n], S[:Vt][1:n, :])
end

"""
    svdvals_restricted(A, Q, n)

Compute the singular values of `A` restricted to the subspace spanned by `Q`
using exact projection.

# Arguments

`A`: input matrix. Must support postmultiply.

`Q`: matrix containing basis vectors of the subspace whose restriction to is
desired.

# Output

`::Vector`: estimated singular values of `A`.

# References

Algorithm 5.1 of \\cite{Halko2011}
"""
function svdvals_restricted(A, Q, n::Int)
    B=Q'A
    S=svdvals!(B)[1:n]
end

"""
    svdfact_re(A, Q)

Compute the SVD factorization of `A` restricted to the subspace spanned by `Q`
using row extraction.

*Note:* \\cite[Remark 5.2]{Halko2011} recommends input of `Q` of the form `Q=A*Ω`
where `Ω` is a sample computed by `randn(n,l)` or even `srft(l)`.

# Arguments

`A`: input matrix. Must support postmultiply

`Q`: matrix containing basis vectors of the subspace whose restriction to is
desired. Need not be orthogonal or normalized.

# Output

`::SVD`: singular value decomposition.

# See also

A faster but less accurate variant of [`svdfact_restricted`](@ref) which uses the
interpolative decomposition `idfact`.

# References

Algorithm 5.2 of \\cite{Halko2011}
"""
function svdfact_re(A, Q)
    F = idfact(Q)
    X, J = F[:B], F[:P]
    R′, W′= qr(A[J, :])
    Z = X*R′
    S=svdfact(Z)
    SVD(S[:U], S[:S], S[:Vt]*W′)
end

"""
    eigfact_restricted(A, Q)

Compute the spectral (`Eigen`) factorization of `A` restricted to the subspace
spanned by `Q` using row extraction.

# Arguments

`A::Hermitian`: input matrix. Must be `Hermitian` and support pre- and post-multiply.

`Q`: orthonormal matrix containing basis vectors of the subspace whose
restriction to is desired.

# Output

`::Base.LinAlg.Eigen`: eigen factorization.

# References

Algorithm 5.3 of \\cite{Halko2011}
"""
function eigfact_restricted(A::Hermitian, Q)
    B = Q'A*Q
    E = eigfact!(B)
    Eigen(E[:values], Q*E[:vectors])
end

"""
    eigfact_re(A, Q)

Compute the spectral (`Eigen`) factorization of `A` restricted to the subspace
spanned by `Q` using row extraction.

*Note:* \\cite[Remark 5.2]{Halko2011} recommends input of `Q` of the form `Q=A*Ω`
where `Ω` is a sample computed by `randn(n,l)` or even `srft(l)`.

# Arguments

`A::Hermitian`: input matrix. Must be `Hermitian` and support pre- and post-multiply.

`Q`: matrix containing basis vectors of the subspace whose restriction to is
desired. Need not be orthogonal or normalized.

# Output

`::Base.LinAlg.Eigen`: eigen factorization.

# See also

A faster but less accurate variant of `eigfact_restricted()` which uses the
interpolative decomposition `idfact()`.

# References

Algorithm 5.4 of \\cite{Halko2011}
"""
function eigfact_re(A::Hermitian, Q)
    X, J = idfact(Q)
    F = qrfact!(X)
    V, R = F[:Q], F[:R]
    Z=R*A[J, J]*R'
    E=eigfact(Z)
    Eigen(E[:values], V*E[:vectors])
end

"""
    eigfact_nystrom(A, Q)

Compute the spectral (`Eigen`) factorization of `A` restricted to the subspace
spanned by `Q` using the Nyström method.

# Arguments

`A`: input matrix. Must be positive semidefinite.

`Q`: orthonormal matrix containing basis vectors of the subspace whose
restriction to is desired.

# Output

`::Base.LinAlg.Eigen`: eigen factorization.

# See also

More accurate than [`eigfact_restricted`](@ref) but is restricted to matrices
that can be Cholesky decomposed.

# References

Algorithm 5.5 of \\cite{Halko2011}
"""
function eigfact_nystrom(A, Q)
    B₁=A*Q
    B₂=Q'*B₁
    C=cholfact!(B₂)
    F=B₁*inv(C)
    S=svdfact!(F)
    Eigen(S[:S].^2, S[:U])
end

"""
    eigfact_onepass(A, Ω)

Compute the spectral (`Eigen`) factorization of `A` using only one matrix
product involving `A`.

# Arguments

`A::Hermitian`: input matrix.

`Ω`: sample matrix for the column space, e.g. `randn(n, l)` or `srft(l)`.

`Ω̃;`: sample matrix for the row space. Not neeeded for `Hermitian` matrices.

`At = A'`: computes transpose of input matrix.

# Output

`::Base.LingAlg.Eigen`: eigen factorization.

# References

Algorithm 5.6 of \\cite{Halko2011}
"""
function eigfact_onepass(A::Hermitian, Ω)
    Y=A*Ω; Q = full(qrfact!(Y)[:Q])
    B=(Q'Y)\(Q'Ω)
    E=eigfact!(B)
    Eigen(E[:values], Q*E[:vectors])
end
function eigfact_onepass(A, Ω, Ω̃; At=A')
    Y=A *Ω; Q = full(qrfact!(Y)[:Q])
    Ỹ=At*Ω; Q̃ = full(qrfact!(Ỹ)[:Q])
    #Want least-squares solution to (5.14 and 5.15)
    B=(Q'Y)\(Q̃'Ω)
    B̃=(Q̃'Ỹ)\(Q'Ω̃)
    #Here is a very very very hacky way to solve the problem
    B=0.5(B + B̃')
    E=eigfact!(B)
    Eigen(E[:values], Q*E[:vectors])
end

"""
    reig(A, l)

Compute the spectral (`Eigen`) decomposition of `A` using a randomized
algorithm.

# Arguments

`A`: input matrix.

`l::Int`: number of eigenpairs to find.

# Output

`::Base.LinAlg.Eigen`: eigen decomposition.

# Implementation note

This is a wrapper around `eigfact_onepass()` which uses the randomized
samples found using `srft(l)`.

# References

```bibtex
@article{Halko2011,
    author = {Halko, N and Martinsson, P G and Tropp, J A},
    doi = {10.1137/090771806},
    journal = {SIAM Review},
    month = jan,
    number = {2},
    pages = {217--288},
    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},
    volume = {53},
    year = {2011}
}
```
"""
reig(A::Hermitian, l::Int) = eigfact_onepass(A, srft(l))
reig(A, l::Int) = eigfact_onepass(A, srft(l), srft(l))
