#############################################
## Randomized singular value decomposition ##
#############################################

# This file provides a rudimentary implementation of the randomized singular
# value decomposition as described in the Reference.
# 
# Reference:
#@article{Halko2011,
#    author = {Halko, N and Martinsson, P G and Tropp, J A},
#    doi = {10.1137/090771806},
#    journal = {SIAM Review},
#    month = jan,
#    number = {2},
#    pages = {217--288},
#    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},
#    volume = {53},
#    year = {2011}
#}


import Base.LinAlg: SVD

export rsvd

"""
Computes the singular value decomposition of "A" using a randomized algorithm.

Inputs:

    A: input matrix
    n: Number of singular value/vector pairs to find
    p: Number of extra vectors to include in computation

Output:

    F: An SVD factorization object

Warning:

    This variant of the randomized singular value decomposition is the most
    commonly found implementation but is not recommended for accurate
    computations, as it often has trouble finding the n largest singular pairs,
    but rather finds n large singular pairs which may not necessarily be the
    largest.

Implementation note:

    This function calls rrange(), which uses naive randomized rangefinding to
    compute a basis for a subspace of dimension n (Algorithm 4.1 of the
    Reference), followed by svdfact_restricted(), which computes the exact SVD
    factorization on the restriction of A to this randomly selected subspace
    (Algorithm 5.1).
"""
function rsvd(A, n, p=0)
    Q = rrange(A, n, p=p)
    svdfact_restricted(A, Q)
end

"""
Computes a basis for a subspace of A of dimension l using naive randomized
rangefinding.

Inputs:

    A: Input matrix. Must support size(A) and premultiply
    l: The number of basis vectors to compute
    p: The oversampling parameter. The number of extra basis vectors to use
       in the computation, which get discarded at the end.

Output:
    Q: A dense matrix of dimension size(A,1) x l containing the basis vectors
       of the computed subspace of A

Warning:

    The Reference explicitly discourages using this algorithm.

Implementation note:

    Whereas the Reference recommends classical Gram-Schmidt with double
    reorthogonalization, we instead compute the basis with qrfact(), which
    for dense A computes the QR factorization using Householder reflectors.
"""
function rrange(A, l::Int; p::Int=0)
    p≥0 || error()
    m, n = size(A)
    if l > m
	warn("Cannot find $l linearly independent vectors of $m x $n matrix")
	warn("Truncating to l=$m, p=0")
        l=m
	p=0
    end
    Ω = randn(n, l+p)
    Y = A*Ω
    Q = full(qrfact!(Y)[:Q])
    Q = p==0 ? Q : Q[:,1:l]
    @assert l==size(Q, 2)
    Q
end


"""
Computes a basis for a subspace of A of dimension l using randomized rangefinding
by subspace iteration.

Inputs:

    A: Input matrix. Must support size(A) and premultiply
    l: The number of basis vectors to compute
    At: (Optional) the transpose of A. Defaults to A'.
    p: The oversampling parameter. The number of extra basis vectors to use
       in the computation, which get discarded at the end. Default: 0
    q: Number of subspace iterations. (Default: 0, which reduces to rrange())

Output:
    Q: A dense matrix of dimension size(A,1) x l containing the basis vectors
       of the computed subspace of A

Implementation note:

    Whereas the Reference recommends classical Gram-Schmidt with double
    reorthogonalization, we instead compute the basis with qrfact(), which
    for dense A computes the QR factorization using Householder reflectors.
"""
function rrange_si(A, l::Integer; At=A', p::Int=0, q::Int=0)
    const basis=_->full(qrfact(_)[:Q])
    p>=0 || error()
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
    p==0 ? Q : Q[:,1:l]
end

"""
Computes the exact SVD factorization of A when restricted to the subspace
spanned by Q.

Inputs:
    A: Input matrix. Must support postmultiply
    Q: Matrix containing basis vectors of the subspace whose restriction to is
       desired

Output:
    F: An SVD Factorization object
"""
function svdfact_restricted(A, Q)
    B=Q'A
    S=svdfact!(B)
    SVD(Q*S[:U], S[:S], S[:Vt])
end

