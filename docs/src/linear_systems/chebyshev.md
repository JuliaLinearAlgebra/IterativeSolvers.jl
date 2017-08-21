# [Chebyshev iteration](@id Chebyshev)

Chebyshev iteration solves the problem $Ax=b$ approximately for $x$ where $A$ is a symmetric, definite linear operator and $b$ the right-hand side vector. The methods assumes the interval $[\lambda_{min}, \lambda_{max}]$ containing all eigenvalues of $A$ is known, so that $x$ can be iteratively constructed via a Chebyshev polynomial with zeros in this interval. This polynomial ultimately acts as a filter that removes components in the direction of the eigenvectors from the initial residual.

The main advantage with respect to Conjugate Gradients is that BLAS1 operations such as inner products are avoided.

## Usage

```@docs
chebyshev
chebyshev!
```

## Implementation details

!!! warning "BLAS1 operations"
    Although the method is often used to avoid computation of inner products, the stopping criterion is still based on the residual norm. Hence the current implementation is not free of BLAS1 operations.

!!! tip
    Chebyshev iteration can be used as an [iterator](@ref Iterators).