# IterativeSolvers.jl

IterativeSolvers.jl is a Julia package that provides iterative algorithms for
solving linear systems, eigenproblems, and singular value problems. The purpose
of this package is to provide efficient Julia implementations for iterative
methods. The package aims to accept a wide variety of input types and that's
why most arguments don't specify a specific type.

For bug reports, feature requests and questions please submit an issue.
If you're interested in contributing, please see the [Contributing](@ref) guide.

For more information on future methods have a look at the package [roadmap](https://github.com/JuliaLang/IterativeSolvers.jl/issues/1) on deterministic methods, for randomized algorithms check [here](https://github.com/JuliaLang/IterativeSolvers.jl/issues/33).

## Linear Solvers

**Stationary methods**

* Jacobi
* Gauss-Seidel
* Successive over-relaxation (SOR)
* Symmetric successive over-relaxation (SSOR)

**Non stationary methods**

* Conjugate Gradients (CG)
* MINRES
* BiCGStabl(l)
* IDR(s)
* Restarted GMRES
* Chebyshev iteration
* LSMR
* LSQR

## Eigenproblem Solvers

* (Inverse) power iteration

## Singular Value Decomposition

* Golub-Kahan-Lanczos
* Randomized singular value decomposition

## Randomized

* Condition number estimate
* Extremal eigenvalue estimates
* Norm estimate

## Documentation Outline

### Manual

```@contents
Pages = [
    "user_manual.md",
]
Depth = 2
```

### Library

```@contents
Pages = ["library/public.md", "library/internal.md"]
Depth = 2
```

### [Index](@id main-index)

### Functions

```@index
Pages = ["library/public.md", "library/internals.md"]
Order = [:function]
```

### Types

```@index
Pages = ["library/public.md", "library/internals.md"]
Order = [:type]
```
