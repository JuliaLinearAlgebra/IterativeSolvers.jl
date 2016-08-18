#########################################################
## Randomized singular value and eigen- decompositions ##
########################################################

# This file provides a rudimentary implementation of the randomized singular
# value decomposition and spectral (eigen-) decomposition as described in
# \cite{Halko2011}.
#
# Reference:
# @article{Halko2011,
#    author = {Halko, N and Martinsson, P G and Tropp, J A},
#    doi = {10.1137/090771806},
#    journal = {SIAM Review},
#    month = jan,
#    number = {2},
#    pages = {217--288},
#    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},
#    volume = {53},
#    year = {2011}
# }

import Base.LinAlg: Eigen, SVD

export rsvdfact, reig

@doc doc"""
Computes the partial singular value decomposition of `A` using a randomized
algorithm.

Inputs:

    `A`: input matrix
    `n`: Number of singular value/vector pairs to find
    `p`: Number of extra vectors to include in computation

Output:

    `F`: An `SVD` `Factorization` object

Warning:

    This variant of the randomized singular value decomposition is the most
    commonly found implementation but is not recommended for accurate
    computations, as it often has trouble finding the `n` largest singular pairs,
    but rather finds `n` large singular pairs which may not necessarily be the
    largest.

Implementation note:

    This function calls `rrange()`, which uses naive randomized rangefinding to
    compute a basis for a subspace of dimension `n` (Algorithm 4.1 of
    \cite{Halko2011}), followed by `svdfact_restricted()`, which computes the
    exact SVD factorization on the restriction of `A` to this randomly selected
    subspace (Algorithm 5.1 of \cite{Halko2011}).

    Alternatively, you can mix and match your own randomized algorithm using
    any of the randomized range finding algorithms to find a suitable subspace
    and feeding the result to one of the routines that computes the SVD
    restricted to that subspace.
""" ->
function rsvdfact(A, n::Int, p::Int=0)
    Q = rrange(A, n+p)
    svdfact_restricted(A, Q, n)
end

@doc doc"""
Like `rsvdfact`, but returns only the singular values.

Inputs:

    as for `rsvdfact`.

Output:

    A vector containing the estimated singular values of `A`

""" ->
function rsvdvals(A, n::Int, p::Int=0)
    Q = rrange(A, n+p)
    svdvals_restricted(A, Q, n)
end

@doc doc"""
Computes an orthonormal basis for a subspace of `A` of dimension `l` using
naive randomized rangefinding.

Inputs:

    `A`: Input matrix. Must support `size(A)` and premultiply
    `l`: The number of basis vectors to compute
    `p`: The oversampling parameter. The number of extra basis vectors to use
       in the computation, which get discarded at the end.

Output:

    `Q`: A dense matrix of dimension `size(A,1)` x l containing the basis
    vectors of the computed subspace of `A`

Reference:

    Algorithm 4.1 of \cite{Halko2011}

Warning:

    The Reference explicitly discourages using this algorithm.

Implementation note:

    Whereas \cite{Halko2011} recommends classical Gram-Schmidt with double
    reorthogonalization, we instead compute the basis with qrfact(), which
    for dense A computes the QR factorization using Householder reflectors.
""" ->
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

@doc doc"""
Computes an orthogonal basis for a subspace of `A` of dimension `l` using adaptive
randomized rangefinding.

Similar to `rrange()`, but determines the oversampling parameter adaptively given
a threshold `ϵ`.

Inputs:

    `A`: Input matrix. Must support `size(A)` and premultiply
    `l`: The number of basis vectors to compute
    `ϵ`: Threshold to determine adaptive fitting
    `maxiter`: Maximum number of iterations to run. (Default: 10.)

Output:

    `Q`: A dense matrix of dimension `size(A,1)` x l containing the basis
         vectors of the computed subspace of `A`

Reference:

    Algorithm 4.2 of \cite{Halko2011}

""" ->
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



@doc doc"""
Computes an orthonormal basis for a subspace of `A` of dimension `l` using
randomized rangefinding by subspace iteration.

Inputs:

    `A` : Input matrix. Must support `size(A)` and premultiply
    `l` : The number of basis vectors to compute
    `At`: (Optional) the transpose of `A`. Defaults to `A'`.
    `p` : The oversampling parameter. The number of extra basis vectors to use
          in the computation, which get discarded at the end. (Default: 0)
    `q` : Number of subspace iterations. (Default: 0)

Note:
    Running with `q=0` is functionally equivalent to `rrange()`

Output:
    `Q`: A dense matrix of dimension `size(A,1)` x `l` containing the basis
         vectors of the computed subspace of `A`

Reference:

    Algorithm 4.4 of \cite{Halko2011}

Implementation note:

    Whereas the Reference recommends classical Gram-Schmidt with double
    reorthogonalization, we instead compute the basis with `qrfact()`, which
    for dense A computes the QR factorization using Householder reflectors.
""" ->
function rrange_si(A, l::Int; At=A', q::Int=0)
    const basis=_->full(qrfact(_)[:Q])
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



@doc doc"""
Computes an orthonormal basis for a subspace of `A` of dimension `l` using
naive randomized rangefinding using stochastic randomized Fourier transforms.

Inputs:

    `A` : Input matrix. Must support `size(A)` and premultiply
    `l` : The number of basis vectors to compute
    `p` : The oversampling parameter. The number of extra basis vectors to use
          in the computation, which get discarded at the end. (Default: 0)

Output:
    `Q`: A dense matrix of dimension `size(A,1)` x `l` containing the basis
         vectors of the computed subspace of `A`

Note:

    Similar to `rrange()`, but does not use Gaussian random matrices.

Reference:

    Algorithm 4.5 of \cite{Halko2011}

Implementation note:

    Whereas the Reference recommends classical Gram-Schmidt with double
    reorthogonalization, we instead compute the basis with `qrfact()`, which
    for dense A computes the QR factorization using Householder reflectors.

""" ->
function rrange_f(A, l::Int)
    n = size(A, 2)
    Ω = srft(l+p)
    Y = A*Ω
    Q = full(qrfact!(Y)[:Q])
end



@doc doc"""
Computes the SVD factorization of `A` restricted to the subspace spanned by `Q`
using exact projection.

Inputs:

    `A`: Input matrix. Must support postmultiply
    `Q`: Matrix containing basis vectors of the subspace whose restriction to is
       desired

Output:

    `F`: An `SVD` `Factorization` object

Reference:

    Algorithm 5.1 of \cite{Halko2011}
""" ->
function svdfact_restricted(A, Q, n::Int)
    B=Q'A
    S=svdfact!(B)
    SVD((Q*S[:U])[:, 1:n], S[:S][1:n], S[:Vt][1:n, :])
end

@doc doc"""
Like `svdfact_restricted`, but returns only the singular values.

Inputs:

    as for `svdfact_restricted`.

Output:

    `v`: A vector containing the estimated singular values of `A`
""" ->
function svdvals_restricted(A, Q, n::Int)
    B=Q'A
    S=svdvals!(B)[1:n]
end

@doc doc"""
Computes the SVD factorization of `A` restricted to the subspace spanned by `Q`
using row extraction.

Inputs:

    `A`: Input matrix. Must support postmultiply
    `Q`: Matrix containing basis vectors of the subspace whose restriction to is
         desired. Need not be orthogonal or normalized.

Note:

    \cite[Remark 5.2]{Halko2011} recommends input of `Q` of the form `Q=A*Ω`
    where `Ω` is a sample computed by `randn(n,l)` or even `srft(l)`.

Output:

    `F`: An `SVD` `Factorization` object

Note:

    A faster but less accurate variant of `svdfact_restricted()` which uses the
    interpolative decomposition `idfact()`.

Reference:

    Algorithm 5.2 of \cite{Halko2011}
""" ->
function svdfact_re(A, Q)
    F = idfact(Q)
    X, J = F[:B], F[:P]
    R′, W′= qr(A[J, :])
    Z = X*R′
    S=svdfact(Z)
    SVD(S[:U], S[:S], S[:Vt]*W′)
end

@doc doc"""
Computes the spectral (`Eigen`) factorization of `A` restricted to the subspace
spanned by `Q` using row extraction.

Inputs:

    `A`: Input matrix. Must be `Hermitian` and support pre- and post-multiply
    `Q`: Orthonormal matrix containing basis vectors of the subspace whose
         restriction to is desired.

Output:

    `F`: An `Eigen` `Factorization` object

Reference:

    Algorithm 5.3 of \cite{Halko2011}
""" ->
function eigfact_restricted(A::Hermitian, Q)
    B = Q'A*Q
    E = eigfact!(B)
    Eigen(E[:values], Q*E[:vectors])
end

@doc doc"""
Computes the spectral (`Eigen`) factorization of `A` restricted to the subspace
spanned by `Q` using row extraction.

Inputs:

    `A`: Input matrix. Must be `Hermitian` and support pre- and post-multiply
    `Q`: Matrix containing basis vectors of the subspace whose restriction to is
         desired. Need not be orthogonal or normalized.

Note:

    \cite[Remark 5.2]{Halko2011} recommends input of `Q` of the form `Q=A*Ω`
    where `Ω` is a sample computed by `randn(n,l)` or even `srft(l)`.

Output:

    `F`: An `Eigen` `Factorization` object

Note:

    A faster but less accurate variant of `eigfact_restricted()` which uses the
    interpolative decomposition `idfact()`.

Reference:

    Algorithm 5.4 of \cite{Halko2011}
""" ->
function eigfact_re(A::Hermitian, Q)
    X, J = idfact(Q)
    F = qrfact!(X)
    V, R = F[:Q], F[:R]
    Z=R*A[J, J]*R'
    E=eigfact(Z)
    Eigen(E[:values], V*E[:vectors])
end

@doc doc"""
Computes the spectral (`Eigen`) factorization of `A` restricted to the subspace
spanned by `Q` using the Nyström method.

Inputs:

    `A`: Input matrix. Must be positive semidefinite.
    `Q`: Orthonormal matrix containing basis vectors of the subspace whose
         restriction to is desired.

Output:

    `F`: An `Eigen` `Factorization` object

Note:

    More accurate than `eigfact_restricted()` but is restricted to matrices
    that can be Cholesky decomposed.

Reference:

    Algorithm 5.5 of \cite{Halko2011}
""" ->
function eigfact_nystrom(A, Q)
    B₁=A*Q
    B₂=Q'*B₁
    C=cholfact!(B₂)
    F=B₁*inv(C)
    S=svdfact!(F)
    Eigen(S[:S].^2, S[:U])
end

@doc doc"""
Computes the spectral (`Eigen`) factorization of `A` using only one matrix
product involving `A`.

Inputs:

    `A`: Input matrix.
    `Ω`: Sample matrix for the column space, e.g. `randn(n, l)` or `srft(l)`
    `Ω̃;: Sample matrix for the row space. Not neeeded for `Hermitian` matrices
    `At`: Computes transpose of input matrix. Default: `A'`

Output:

    `F`: An `Eigen` `Factorization` object

Reference:

    Algorithm 5.6 of \cite{Halko2011}
""" ->
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

@doc doc"""
Computes the spectral (`Eigen`) decomposition of `A` using a randomized
algorithm.

Inputs:

    `A`: input matrix
    `n`: Number of eigenpairs to find

Output:

    `F`: An `Eigen` `Factorization` object

Implementation Note:

    This is a wrapper around `eigfact_onepass()` which uses the randomized
    samples found using `srft(l)`.
""" ->
reig(A::Hermitian, l::Int) = eigfact_onepass(A, srft(l))
reig(A, l::Int) = eigfact_onepass(A, srft(l), srft(l))
