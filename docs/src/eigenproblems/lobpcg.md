# [Locally optimal block preconditioned conjugate gradient (LOBPCG)](@id LOBPCG)

Solves the generalized eigenproblem $Ax = λBx$ approximately where $A$ and $B$ are Hermitian linear maps, and $B$ is positive definite. $B$ is taken to be the identity by default. It can find the smallest (or largest) `k` eigenvalues and their corresponding eigenvectors which are B-orthonormal. It also admits a preconditioner and a "constraints" matrix `C`, such that the algorithm returns the smallest (or largest) eigenvalues associated with the eigenvectors in the nullspace of `C'B`.

## Usage

```@docs
lobpcg
lobpcg!
```

## Implementation Details

A `LOBPCGIterator` is created to pre-allocate all the memory required by the method using the constructor `LOBPCGIterator(A, B, largest, X, P, C)` where `A` and `B` are the matrices from the generalized eigenvalue problem, `largest` indicates if the problem is a maximum or minimum eigenvalue problem, `X` is the initial eigenbasis, randomly sampled if not input, where `size(X, 2)` is the block size `bs`. `P` is the preconditioner, `nothing` by default, and `C` is the constraints matrix. The desired `k` eigenvalues are found `bs` at a time.


## References

[Knyazev1993]:

    Andrew V. Knyazev. "Toward the Optimal Preconditioned Eigensolver: Locally
    Optimal Block Preconditioned Conjugate Gradient Method" SIAM Journal
    on Scientific Computing, 23(2):517–541, 2001.

[Scipy LOBPCG implementation]:

    https://github.com/scipy/scipy/blob/v1.1.0/scipy/sparse/linalg/eigen/lobpcg/lobpcg.py#L109-L568
