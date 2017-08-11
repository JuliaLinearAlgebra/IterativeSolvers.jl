# [MINRES](@id MINRES)

MINRES is a short-recurrence version of GMRES for solving $Ax = b$ approximately for $x$ where $A$ is a symmetric, Hermitian, skew-symmetric or skew-Hermitian linear operator and $b$ the right-hand side vector.

## Usage

```@docs
minres
minres!
```

## Implementation details
MINRES exploits the tridiagonal structure of the Hessenberg matrix. Although MINRES is mathematically equivalent to GMRES, it might not be equivalent in finite precision. MINRES updates the solution as

$$x := x_0 + (V R^{-1}) (Q^*\|r_0\|e_1)$$

where $V$ is the orthonormal basis for the Krylov subspace and $QR$ is the QR-decomposition of the Hessenberg matrix. Note that the brackets are placed slightly differently from how GMRES would update the residual.

MINRES computes $V$ and $W = VR^{-1}$ via a three-term recurrence, using only the last column of $R.$ Therefore we pre-allocate only six vectors, save only the last two entries of $Q^*\|r_0\|e_1$ and part of the last column of the Hessenberg matrix.

!!! note "Real and complex arithmetic"
    If $A$ is Hermitian, then the Hessenberg matrix will be real. This is exploited in the current implementation.

    If $A$ is skew-Hermitian, the diagonal of the Hessenberg matrix will be imaginary, and hence we use complex arithmetic in that case.
