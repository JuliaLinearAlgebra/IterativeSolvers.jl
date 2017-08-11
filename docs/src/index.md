# IterativeSolvers.jl

IterativeSolvers.jl is a Julia package that provides efficient iterative algorithms for solving large linear systems, eigenproblems, and singular value problems. Most of the methods can be used *matrix-free*.

For bug reports, feature requests and questions please submit an issue. If you're interested in contributing, please see the [Contributing](@ref) guide.

For more information on future methods have a look at the package [roadmap](https://github.com/JuliaLang/IterativeSolvers.jl/issues/1) on deterministic methods, for randomized algorithms check [here](https://github.com/JuliaLang/IterativeSolvers.jl/issues/33).

## What method should I use for linear systems?

When solving linear systems $Ax = b$ for a square matrix $A$ there are quite some options. The typical choices are listed below:

| Method              | When to use it                                                           |
|---------------------|--------------------------------------------------------------------------|
| [Conjugate Gradients](@ref CG) | Best choice for **symmetric**, **positive-definite** matrices |
| [MINRES](@ref MINRES) | For **symmetric**, **indefinite** matrices |
| [GMRES](@ref GMRES) | For **nonsymmetric** matrices when a good preconditioner is available |
| [IDR(s)](@ref IDRs) | For **nonsymmetric**, **strongly indefinite** problems without a good preconditioner |
| [BiCGStab(l)](@ref BiCGStabl) | Otherwise for **nonsymmetric** problems |

We also offer [Chebyshev iteration](@ref Chebyshev) as an alternative to Conjugate Gradients when bounds on the spectrum are known.

Stationary methods like [Jacobi](@ref), [Gauss-Seidel](@ref), [SOR](@ref) and [SSOR](@ref) can be used as smoothers to reduce high-frequency components in the error in just a few iterations.

When solving **least-squares** problems we currently offer just [LSMR](@ref LSMR) and [LSQR](@ref LSQR).

## Eigenproblems and SVD

For the Singular Value Decomposition we offer [SVDL](@ref SVDL), which is the Golub-Kahan-Lanczos procedure.

For eigenvalue problems we have at this point just the [Power Method](@ref PowerMethod) and some convenience wrappers to do shift-and-invert.

## Randomized algorithms

[Randomized algorithms](@ref Randomized) have gotten some traction lately. Some of the methods mentioned in [^Halko2011] have been implemented as well, although their quality is generally poor. Also note that many classical methods such as the subspace iteration, BiCG and recent methods like IDR(s) are also "randomized" in some sense.

[^Halko2011]: Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions." SIAM review 53.2 (2011): 217-288.