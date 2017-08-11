# [Golub-Kahan-Lanczos (SVDL)](@id SVDL)

The SVDL method computes a partial, approximate SVD decomposition of a general linear operator $A$.

## Usage

```@docs
svdl
```

## Implementation details

The implementation of thick restarting follows closely that of SLEPc as described in [^Hernandez2008]. Thick restarting can be turned off by setting `k = maxiter`, but most of the time this is not desirable.

The singular vectors are computed directly by forming the Ritz vectors from the product of the Lanczos vectors `L.P`/`L.Q` and the singular vectors of `L.B`. Additional accuracy in the singular triples can be obtained using inverse iteration.

[^Hernandez2008]: Vicente Hernández, José E. Román, and Andrés Tomás. "A robust and efficient parallel SVD solver based on restarted Lanczos bidiagonalization." Electronic Transactions on Numerical Analysis 31 (2008): 68-85.
