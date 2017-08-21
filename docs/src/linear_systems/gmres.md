# [Restarted GMRES](@id GMRES)

GMRES solves the problem $Ax = b$ approximately for $x$ where $A$ is a general, linear operator and $b$ the right-hand side vector. The method is optimal in the sense that it selects the solution with minimal residual from a Krylov subspace, but the price of optimality is increasing storage and computation effort per iteration. Restarts are necessary to fix these costs.

## Usage

```@docs
gmres
gmres!
```

## Implementation details

The implementation pre-allocates a matrix $V$ of size `n` by `restart` whose columns form an orthonormal basis for the Krylov subspace. This allows BLAS2 operations when updating the solution vector $x$. The Hessenberg matrix is also pre-allocated.

Modified Gram-Schmidt is used to orthogonalize the columns of $V$.

The computation of the residual norm is implemented in a non-standard way, namely keeping track of a vector $\gamma$ in the null-space of $H_k^*$, which is the adjoint of the $(k + 1) \times k$ Hessenberg matrix $H_k$ at the $k$th iteration. Only when $x$ needs to be updated is the Hessenberg matrix mutated with Givens rotations.

!!! tip
    GMRES can be used as an [iterator](@ref Iterators). This makes it possible to access the Hessenberg matrix and Krylov basis vectors during the iterations.