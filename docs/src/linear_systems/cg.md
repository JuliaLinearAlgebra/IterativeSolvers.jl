# [Conjugate Gradients (CG)](@id CG)

Conjugate Gradients solves $Ax = b$ approximately for $x$ where $A$ is a symmetric, positive-definite linear operator and $b$ the right-hand side vector. The method uses short recurrences and therefore has fixed memory costs and fixed computational costs per iteration.

## Usage

```@docs
cg
cg!
```

## Implementation details

The current implementation follows a rather standard approach. Note that preconditioned CG (or PCG) is slightly different from ordinary CG, because the former must compute the residual explicitly, while it is available as byproduct in the latter. Our implementation of CG ensures the minimal number of vector operations.

!!! tip
    CG can be used as an [iterator](@ref Iterators).