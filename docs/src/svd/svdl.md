# [Golub-Kahan-Lanczos (SVDL)](@id SVDL)

The SVDL method computes a partial, approximate SVD decomposition of a general linear operator $A$.

## Usage

```@docs
svdl
```

## Implementation details

The implementation of thick restarting follows closely that of SLEPc as described in [^Hernandez2008]. Thick restarting can be turned off by setting `k = maxiter`, but most of the time this is not desirable.

The singular vectors are computed directly by forming the Ritz vectors from the product of the Lanczos vectors `L.P`/`L.Q` and the singular vectors of `L.B`. Additional accuracy in the singular triples can be obtained using inverse iteration.

A deterministic seed is used for generating pseudo-random initial
data for the algorithm; this can be controlled by passing a
different pseudorandom number generator (an [`AbstractRNG`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG)) via
the `rng` keyword argument, or by passing an initial `v0` vector
directly.

[^Hernandez2008]: Vicente Hernández, José E. Román, and Andrés Tomás. "A robust and efficient parallel SVD solver based on restarted Lanczos bidiagonalization." Electronic Transactions on Numerical Analysis 31 (2008): 68-85.
