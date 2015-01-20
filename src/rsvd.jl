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

@doc doc"""
Computes the partial singular value decomposition of `A` using a randomized
algorithm.

Inputs:

    `A`: input matrix
    `n`: Number of singular value/vector pairs to find
    `p`: Number of extra vectors to include in computation

Output:

    `F`: An SVD factorization object

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
""" ->
function rsvd(A, n, p=0)
    Q = rrange(A, n, p=p)
    svdfact_restricted(A, Q)
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

@doc doc"""
Computes an orthogonal basis for a subspace of A of dimension l using adaptive
randomized rangefinding.

Similar to rrange(), but determines the oversampling parameter adaptively given
a threshold ϵ.

Inputs:

    `A`: Input matrix. Must support `size(A)` and premultiply
    `l`: The number of basis vectors to compute
    `ϵ`: Threshold to determine adaptive fitting
    `maxiter`: Maximum number of iterations to run. Default: 10.

Output:

    `Q`: A dense matrix of dimension `size(A,1)` x l containing the basis
         vectors of the computed subspace of `A`

Reference:

    Algorithm 4.2 of \cite{Halko2011}

Implementation note:

    Whereas \cite{Halko2011} recommends classical Gram-Schmidt with double
    reorthogonalization, we instead compute the basis with `qrfact()`, which
    for dense `A` computes the QR factorization using Householder reflectors.

Warning:

    Empirical testing indicates that this implementation is slow; the
    orthogonalization of a newly computed random vector is rate-limiting.
""" ->
function rrange_adaptive(A, r::Integer, ϵ::Real=eps(); maxiter::Int=10)
    m, n = size(A)
    r += maxiter
    Ω = randn(n, r)
    #Normalize columns of Ω
    #for i=1:r 
    #    ω = sub(Ω, :, i)
    #    scale!(1/norm(ω), ω)
    #end
    Y = A*Ω
    Q = full(qrfact(Y)[:Q], thin=true)

    const tol=ϵ/(10*√(2/π))
    for j=1:maxiter
        Tol = maximum([norm(sub(Y,:,i)) for i=j-1+(1:r)])
        info("Iteration $j: norm = $Tol (target: $tol)")
        Tol > tol || break

        if j>1
            y = sub(Y, :, j)
            Y[:,j]-=Q*(Q'y)
        end
        q = sub(Y,:,j)
        scale!(1/norm(q), q)
        Q = [Q q]

	#Compute new random vector in the range of A and
	#orthogonalize against existing vectors
        y = randn(n)
        scale!(1/norm(y), y)
        A_mul_B!(y, A, y)
        y -= Q*Q'y
        Y = [Y y]
        for i = j+(1:(r-1))
            y = sub(Y, :, i)
            Y[:,i] -= q*(q⋅y)
        end
        j==maxiter && warn("Maximum number of iterations reached")
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

@doc doc"""
Computes the exact SVD factorization of `A` when restricted to the subspace
spanned by `Q`.

Inputs:

    `A`: Input matrix. Must support postmultiply
    `Q`: Matrix containing basis vectors of the subspace whose restriction to is
       desired

Output:

    `F`: An `SVD` `Factorization` object

Reference:

    Algorithm 5.1 of \cite{Halko2011}
""" ->
function svdfact_restricted(A, Q)
    B=Q'A
    S=svdfact!(B)
    SVD(Q*S[:U], S[:S], S[:Vt])
end

@doc doc"""
Computes the SVD factorization of `A` when restricted to the subspace
spanned by `Q` using row extraction.

Inputs:

    `A`: Input matrix. Must support postmultiply
    `Q`: Matrix containing basis vectors of the subspace whose restriction to is
         desired. Need not be orthogonal or normalized.

Note:

    \cite{Halko2011} recommends input of `Q` of the form `Q=A*randn(n,l)`.

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
