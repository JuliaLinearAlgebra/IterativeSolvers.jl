# [Golub-Kahan-Lanczos (SVDL)](@id SVDL)

The SVDL method computes a partial, approximate SVD decomposition of a general linear operator $A$.

## Usage

```@docs
svdl
```

## Implementation details

The implementation of thick restarting follows closely that of SLEPc as
described in [Hernandez2008]. Thick restarting can be turned off by setting `k
= maxiter`, but most of the time this is not desirable.

The singular vectors are computed directly by forming the Ritz vectors from the
product of the Lanczos vectors `L.P`/`L.Q` and the singular vectors of `L.B`.
Additional accuracy in the singular triples can be obtained using inverse
iteration.

## References

```bibtex
@article{Golub1965,
    author = {Golub, G. and Kahan, W.},
    doi = {10.1137/0702016},
    journal = {Journal of the Society for Industrial and Applied Mathematics
        Series B Numerical Analysis},
    volume = 2,
    number = 2,
    pages = {205--224},
    title = {Calculating the Singular Values and Pseudo-Inverse of a Matrix},
    year = 1965
}

@article{Wu2000,
    author = {Wu, Kesheng and Simon, Horst},
    journal = {SIAM Journal on Matrix Analysis and Applications},
    number = 2,
    pages = {602--616},
    title = {Thick-Restart {L}anczos Method for Large Symmetric Eigenvalue Problems},
    volume = 22,
    year = 2000
}

@article{Baglama2005,
    author = {Baglama, James and Reichel, Lothar},
    doi = {10.1137/04060593X},
    journal = {SIAM Journal on Scientific Computing},
    number = 1,
    pages = {19--42},
    title = {Augmented Implicitly Restarted {L}anczos Bidiagonalization Methods},
    volume = 27,
    year = 2005
}

@article{Hernandez2008,
    author = {Hern\'{a}ndez, Vicente and Rom\'{a}n, Jos\'{e} E and Tom\'{a}s,
    Andr\'{e}s},
    journal = {Electronic Transactions on Numerical Analysis},
    pages = {68--85},
    title = {A Robust and Efficient Parallel {SVD} Solver based on Restarted
        {L}anczos Bidiagonalization},
    url = {http://etna.mcs.kent.edu/volumes/2001-2010/vol31/abstract.php?vol=31\&pages=68-85},
    volume = 31,
    year = 2008
}
```