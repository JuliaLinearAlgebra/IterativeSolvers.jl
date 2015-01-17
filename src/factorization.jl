############################
# Specialized factorizations
############################

export idfact

immutable Interpolative{T} <: Factorization{T}
    B :: AbstractMatrix{T}
    P :: AbstractMatrix{T}
end

"""
Compute the interpolative decomposition of A

    A ≈ B * P

    B's columns are a subset of the columns of A
    some subset of P's columns are the k x k identity,
    no entry of P exceeds magnitude 2, and
    ||B * P - A|| ≲ σ(A, k+1)

Inputs:
    A: Matrix to factorize
    k: Number of columns of A to return in B
    l: Length of random vectors to project onto

Reference:

    Algorithm I of Liberty2007

Implementation note:

    This is a hacky version of the algorithms described in Liberty2007 and
    Cheng2005. The former refers to the factorization (3.1) of the latter.
    However, it is not actually necessary to compute this factorization in
    its entirely to compute an interpolative decomposition. Instead, it
    suffices to find some permutation of the first k columns of Y = R * A,
    extract the subset of A into B, then compute the P matrix as B\A which
    will automatically compute P using a suitable least-squares algorithm.
    
    The approximation we use here is to compute the column pivots of Y,
    rather then use the true column pivots as would be computed by a column-
    pivoted QR process.

    @article{Cheng2005,
        author = {Cheng, H and Gimbutas, Z and Martinsson, P G and Rokhlin, V},
        doi = {10.1137/030602678},
        issn = {1064-8275},
        journal = {SIAM Journal on Scientific Computing},
        month = jan,
        number = {4},
        pages = {1389--1404},
        title = {On the Compression of Low Rank Matrices},
        volume = {26},
        year = {2005}
    }
"""
function idfact(A, k::Int, l::Int)
    m, n = size(A)
    R = randn(l, m)
    Y = R * A #size l x n

    #Compute column pivots of first k columns of Y
    maxvals = map(j->maximum(abs(sub(Y, :, j))), 1:n)
    piv = sortperm(maxvals, rev=true)[1:k]

    B = A[:, piv]
    Interpolative(B, B\A)
end

