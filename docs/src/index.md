# IterativeSolvers.jl

IterativeSolvers.jl is a Julia package that provides efficient iterative algorithms for solving large linear systems, eigenproblems, and singular value problems. Most of the methods can be used *matrix-free*.

For bug reports, feature requests and questions please submit an issue. If you're interested in contributing, please see the [Contributing](@ref) guide.

For more information on future methods have a look at the package [roadmap](https://github.com/JuliaLang/IterativeSolvers.jl/issues/1).

## What method should I use for linear systems?

### Square linear systems

When solving linear systems $Ax = b$ for a square matrix $A$ there are quite some options. The typical choices are listed below:

| Method              | When to use it                                                           |
|---------------------|--------------------------------------------------------------------------|
| [Conjugate Gradients](@ref CG) | Best choice for **symmetric**, **positive-definite** matrices |
| [MINRES](@ref MINRES) | For **symmetric**, **indefinite** matrices |
| [GMRES](@ref GMRES) | For **nonsymmetric** matrices when a good [preconditioner](@ref Preconditioning) is available |
| [IDR(s)](@ref IDRs) | For **nonsymmetric**, **strongly indefinite** problems without a good preconditioner |
| [BiCGStab(l)](@ref BiCGStabl) | Otherwise for **nonsymmetric** problems |

We also offer [Chebyshev iteration](@ref Chebyshev) as an alternative to Conjugate Gradients when bounds on the spectrum are known.

Stationary methods like [Jacobi](@ref), [Gauss-Seidel](@ref), [SOR](@ref) and [SSOR](@ref) can be used as smoothers to reduce high-frequency components in the error in just a few iterations.

### Non-square systems: least squares

When solving **least-squares** problems we currently offer [LSMR](@ref LSMR) and [LSQR](@ref LSQR).
LSMR generally converges more quickly than LSQR.

## Eigenproblems and SVD

For the Singular Value Decomposition we offer [SVDL](@ref SVDL), which is the Golub-Kahan-Lanczos procedure.

For eigenvalue problems we have at this point just the [Power Method](@ref PowerMethod) and some convenience wrappers to do shift-and-invert.
