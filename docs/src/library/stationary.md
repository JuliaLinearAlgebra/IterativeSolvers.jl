# Stationary methods

Stationary methods are typically used as smoothers in multigrid methods, where only very few iterations are applied to get rid of high-frequency components in the error. The implementations of stationary methods have this goal in mind, which means there is no other stopping criterion besides the maximum number of iterations.

!!! note "CSC versus CSR"
    Julia stores matrices column-major. In order to avoid cache misses, the implementations of our stationary methods traverse the matrices column-major. This deviates from classical textbook implementations. Also the (S)SOR methods cannot be computed fully in-place, but require a temporary vector.

    When it comes to `SparseMatrixCSC`, we precompute an integer array of the indices of the diagonal as well to avoid expensive searches in each iteration.

## Jacobi

```@docs
jacobi
jacobi!
```

## Gauss-Seidel

```@docs
gauss_seidel
gauss_seidel!
```

## Successive over-relaxation (SOR)

```@docs
sor
sor!
```

## Symmetric successive over-relaxation (SSOR)

```@docs
ssor
ssor!
```
