############################
# Specialized factorizations
############################

export idfact

"""
An interpolative decomposition.

For a matrix `A`, the interpolative decomposition `F` contains the matrices `B`
and `P` computed by `idfact()`. See the documentation of `idfact()` for details.

# References

\\cite{Cheng2005, Liberty2007}
"""
struct Interpolative{T} <: Factorization{T}
    B :: AbstractMatrix{T}
    P :: AbstractMatrix{T}
end

"""
    idfact(A, k, l)

Compute and return the interpolative decomposition of `A`: A ≈ B * P

Where:
* `B`'s columns are a subset of the columns of `A`
* some subset of `P`'s columns are the `k x k` identity, no entry of `P` exceeds magnitude 2, and
* ||B * P - A|| ≲ σ(A, k+1), the (`k+1`)st singular value of `A`.

# Arguments

`A`: Matrix to factorize

`k::Int`: Number of columns of A to return in B

`l::Int`: Length of random vectors to project onto

# Output

`(::Interpolative)`: interpolative decomposition.

# Implementation note

This is a hacky version of the algorithms described in \\cite{Liberty2007}
and \\cite{Cheng2005}. The former refers to the factorization (3.1) of the
latter.  However, it is not actually necessary to compute this
factorization in its entirely to compute an interpolative decomposition.

Instead, it suffices to find some permutation of the first k columns of Y =
R * A, extract the subset of A into B, then compute the P matrix as B\\A
which will automatically compute P using a suitable least-squares
algorithm.

The approximation we use here is to compute the column pivots of Y,
rather then use the true column pivots as would be computed by a column-
pivoted QR process.

# References

\\cite[Algorithm I]{Liberty2007}

```bibtex
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
```
"""
function idfact(A, k::Int, l::Int)
    m, n = size(A)
    R = randn(l, m)
    Y = R * A #size l x n

    #Compute column pivots of first k columns of Y
    maxvals = map(j->maximum(abs(view(Y, :, j))), 1:n)
    piv = sortperm(maxvals, rev=true)[1:k]

    B = A[:, piv]
    Interpolative(B, B\A)
end
