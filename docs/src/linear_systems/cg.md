# [Conjugate Gradients (CG)](@id CG)

Conjugate Gradients solves $Ax = b$ approximately for $x$ where $A$ is a symmetric, positive-definite linear operator and $b$ the right-hand side vector. The method uses short recurrences and therefore has fixed memory costs and fixed computational costs per iteration.

## Usage

```@docs
cg
cg!
```

## On the GPU

The method should work fine on the GPU. As a minimal working example, consider:

```julia
using LinearAlgebra, CuArrays, IterativeSolvers

n = 100
A = cu(rand(n, n))
A = A + A' + 2*n*I
b = cu(rand(n))
x = cg(A, b)
```

!!! note
    Make sure that all state vectors are stored on the GPU. For instance when calling `cg!(x, A, b)`, one might have an issue when `x` is stored on the GPU, while `b` is stored on the CPU -- IterativeSolvers.jl does not copy the vectors to the same device.


## Implementation details

The current implementation follows a rather standard approach. Note that preconditioned CG (or PCG) is slightly different from ordinary CG, because the former must compute the residual explicitly, while it is available as byproduct in the latter. Our implementation of CG ensures the minimal number of vector operations.

!!! tip
    CG can be used as an [iterator](@ref Iterators).