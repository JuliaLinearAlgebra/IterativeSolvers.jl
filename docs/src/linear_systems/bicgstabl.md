# [BiCGStab(l)](@id BiCGStabl)

BiCGStab(l) solves the problem $Ax = b$ approximately for $x$ where $A$ is a general, linear operator and $b$ the right-hand side vector. The methods combines BiCG with $l$ GMRES iterations, resulting in a short-recurrence iteration. As a result the memory is fixed as well as the computational costs per iteration.

## Usage

```@docs
bicgstabl
bicgstabl!
```

## Implementation details

The method is based on the original article [^Sleijpen1993], but does not implement later improvements. The normal equations arising from the GMRES steps are solved without orthogonalization. Hence the method should only be reliable for relatively small values of $l$.

The `r` and `u` factors are pre-allocated as matrices of size $n \times (l + 1)$, so that BLAS2 methods can be used. Also the random shadow residual is pre-allocated as a vector. Hence the storage costs are approximately $2l + 3$ vectors.

A deterministic seed is used for generating pseudo-random initial
data for the algorithm; this can be controlled by passing a
different pseudorandom number generator (an [`AbstractRNG`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG)) via
the `rng` keyword argument.

!!! tip
    BiCGStabl(l) can be used as an [iterator](@ref Iterators).

[^Sleijpen1993]:

    Sleijpen, Gerard LG, and Diederik R. Fokkema. "BiCGstab(l) for
    linear equations involving unsymmetric matrices with complex spectrum."
    Electronic Transactions on Numerical Analysis 1.11 (1993): 2000.