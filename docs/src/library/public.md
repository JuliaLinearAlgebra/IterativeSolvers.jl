# Public Documentation

Documentation for `IterativeSolvers.jl`'s public interface.

```@contents
Pages = ["public.md"]
Depth = 4
```

## Index

```@index
Pages = ["public.md"]
```

## Types

### `ConvergenceHistory`

```@docs
ConvergenceHistory
```

### `KrylovSubspace`

```@docs
KrylovSubspace
```

### `FuncMat`

```@docs
FuncMat
```

## Linear Solvers

### `Jacobi`

```@docs
jacobi
jacobi!
```

### `Gauss-Seidel`

```@docs
gauss_seidel
gauss_seidel!
```

### `Successive over relaxation`

```@docs
sor
sor!
```

### `Symmetric successive over relaxation`

```@docs
ssor
ssor!
```

### `IDRS`

```@docs
idrs
idrs!
```
**References**

```
[1] IDR(s): a family of simple and fast algorithms for solving large
    nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen
    SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, 2008
[2] Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits
    Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld
    ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011
[3] This file is a translation of the following MATLAB implementation:
    http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m
[4] IDR(s)' webpage http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html
```

### `LSMR`

```@docs
lsmr
lsmr!
```

**References**

Adapted from: [http://web.stanford.edu/group/SOL/software/lsmr/](http://web.stanford.edu/group/SOL/software/lsmr/)

### `LSQR`

```@docs
lsqr
lsqr!
```

**References**

```
Adapted from: http://web.stanford.edu/group/SOL/software/lsqr/

1. C. C. Paige and M. A. Saunders (1982a).
    LSQR: An algorithm for sparse linear equations and sparse least squares,
    ACM TOMS 8(1), 43-71.

2. C. C. Paige and M. A. Saunders (1982b).
    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,
    ACM TOMS 8(2), 195-209.

3. M. A. Saunders (1995).  Solution of sparse rectangular systems using
    LSQR and CRAIG, BIT 35, 588-604.
```


### `Conjugate gradients`

```@docs
cg
cg!
```

### `Chebyshev iteration`

```@docs
chebyshev
chebyshev!
```

### `Generalized minimal residual method (with restarts)`

```@docs
gmres
gmres!
```

## Eigen Solvers

### `Power iteration`

```@docs
powm
```

### `Inverse power iteration`

```@docs
invpowm
```

### `Simple Lanczos`

```@docs
eiglancz
```

### `Golub-Kahan-Lanczos`

```@docs
svdl
```

**Implementation notes**

The implementation of thick restarting follows closely that of SLEPc as
described in [Hernandez2008]. Thick restarting can be turned off by setting `k
= maxiter`, but most of the time this is not desirable.

The singular vectors are computed directly by forming the Ritz vectors from the
product of the Lanczos vectors `L.P`/`L.Q` and the singular vectors of `L.B`.
Additional accuracy in the singular triples can be obtained using inverse
iteration.

**References**

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


## Randomized

### `Condition number estimate`

```@docs
rcond
```

### `Extremal eigenvalue estimates`

```@docs
reigmin
reigmax
```

### `Norm estimate`

```@docs
rnorm
rnorms
```

### `Randomized singular value decomposition`

```@docs
reig
rsvdfact
rsvd_fnkz
```
