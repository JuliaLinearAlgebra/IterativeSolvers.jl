# Iterative Solvers

[![IterativeSolvers](http://pkg.julialang.org/badges/IterativeSolvers_0.4.svg)](http://pkg.julialang.org/?pkg=IterativeSolvers)
[![IterativeSolvers](http://pkg.julialang.org/badges/IterativeSolvers_0.5.svg)](http://pkg.julialang.org/?pkg=IterativeSolvers)
[![Build Status on Linux](https://travis-ci.org/JuliaMath/IterativeSolvers.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/IterativeSolvers.jl)
[![Coveralls](https://coveralls.io/repos/github/JuliaMath/IterativeSolvers.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/IterativeSolvers.jl?branch=master)
[![Codecov](https://codecov.io/gh/JuliaMath/IterativeSolvers.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/IterativeSolvers.jl)

&copy; 2013---2016 [The Contributors](https://github.com/JuliaLang/IterativeSolvers.jl/contributors). Released under the [MIT License](https://github.com/JuliaLang/julia/blob/master/LICENSE.md).

`IterativeSolvers` is a [Julia](http://julialang.org) package that provides iterative algorithms for solving linear systems, eigensystems, and singular value problems.

## Functions provided

- `svdl`: Compute the largest few singular triples using a iterative singular
          value decomposition.

- `rsvdfact`: Compute the largest few singular triples using a naive randomized
              SVD.

## Current status

- [Issue #1](https://github.com/JuliaLang/IterativeSolvers.jl/issues/1) documents the implementation roadmap for iterative solvers.
- The interfaces between the various algorithms are still in flux, but aim to be consistent.

- [Issue #33](https://github.com/JuliaLang/IterativeSolvers.jl/issues/33) documents the implementation roadmap for randomized algorithms.

## Consistent interface

### Nomenclature

- All randomized algorithms start with `r_`.
- All linear solvers do not have a prefix. Routines that accept preconditioners do *not* take a `p` prefix, e.g. `cg` instead of `pcg`.
- All other eigenvalue solvers start with `eigvals_`.
- All other singular value solvers start with `svdvals_`.

### Interface

- All linear-algebraic routines will take as input a linear operator `A` that maps vectors to vectors. `A` is not explicitly typed, but must either be a `KrylovSubspace` or support multiplication `*` or function composition (`apply`) that behave as necessary to produce the correct mapping on the vector space.

- A custom type for `A` may be specified. The following interface is expected to be defined on `A`:
  - `A*v` is define and computes the matrix-vector product on a `v::Vector`.
  - `eltype(A)` is defined and returns the element type implicit in the equivalent matrix representation of `A`.
  - `size(A, d)` is defined and returns the nominal dimensions along the `d`th axis in the equivalent matrix representation of `A`.

- All linear solvers have a common function declaration

    solver(A, b::Vector, [x, Pl, Pr, varargs...]; tol::Real, maxiter::Int, kwargs...)

  - `A` is a `KrylovSubspace` or linear operator as described above.
  - `b` is the vector to be solved
  - `x` is a vector for the initial guess (if not specified, defaults to random unit vector on the unit hypersphere)
  - `Pl` is the (left) preconditioner (for routines that accept them)
  - `Pr` is the right preconditioner (for routines that accept them)
  - `tol` is the threshold for determining convergence (which is typically compared to the norm of the residual vector, but whose semantics may vary from solver to solver) and conventionally defaults to `size(A,1)^3*eps()`
  - `maxiter` is the maximum number of allowed iterations, which conventionally defaults to `length(b)`
  - additional `varargs` and `kwargs` as necessary.

- All linear solvers have a common return format
  - `x::Vector`, the solution vector
  - `h::ConvergenceHistory`, a special data type returning the convergence history as discussed in [Issue #6](https://github.com/JuliaLang/IterativeSolvers.jl/issues/6)

### Krylov subspaces

The `KrylovSubspace` type collects information on the [Krylov subspace](http://en.wikipedia.org/wiki/Krylov_subspace) generated over the course of an iterative Krylov solver.

Recall that the Krylov subspace of order `r` given a starting vector `b` and a linear operator `A` is spanned by the vectors `[b, A*b, A^2*b,... A^(r-1)*b]`. Many modern iterative solvers work on Krylov spaces which expand until they span enough of the range of `A` for the solution vector to converge. Most practical algorithms, however, will truncate the order of the Krylov subspace to save memory and control the accumulation of roundoff errors. Additionally, they do not work directly with the raw iterates `A^n*b`, but will orthogonalize subsequent vectors as they are computed so as to improve numerical stability. `KrylovSubspace`s provide a natural framework to express operations such as a (possibly non-orthonormalized) basis for the Krylov subspace, retrieving the next vector in the subspace, and orthogonalizing an arbitrary vector against (part or all of) the subspace.

The implementation of `KrylovSubspace` in this package differs from standard textbook presentations of iterative solvers. First, the `KrylovSubspace` type shows clearly the relationship between the linear operator `A` and the sequence of basis vectors for the Krylov subspace that is generated by each new iteration. Second, the grouping together of basis vectors also makes the orthogonalization steps in each iteration routine clearer. Third, unlike in other languages, the powerful type system of Julia allows for a simplified notation for iterative algorithms without compromising on performance, and enhances code reuse across multiple solvers.

#### Constructors

A `KrylovSubspace` can be initialized by three constructors, depending on the type of the linear operator:

    KrylovSubspace{T}(A::AbstractMatrix{T}, [order::Int, v::Vector{Vector{T}}])

    KrylovSubspace{T}(A::KrylovSubspace{T}, [order::Int, v::Vector{Vector{T}}])

    KrylovSubspace(A, n::Int, [order::Int, v::Vector{Vector{T}}])

- `A`: The linear operator associated with the `KrylovSubspace`.
- `order`: the order of the `KrylovSubspace`, i.e. the maximal number of Krylov vectors to remember.
- `n`: the dimensionality of the underlying vector space `T^n`.
- `v`: the iterable collection of Krylov vectors (of maximal length `order`).

  The dimensionality of the underlying vector space is automatically inferred where possible, i.e. when the linear operator is an `AbstractMatrix` or `KrylovSubspace`. (Note: the second constructor destroys the old `KrylovSubspace`.)

### Convergence history

A [`ConvergenceHistory`](https://github.com/JuliaLang/IterativeSolvers.jl/issues/6) type provides information about the iteration history. Currently, `ConvergenceHistory` provides:
  - `isconverged::Bool`, a flag for whether or not the algorithm is converged
  - `threshold`, the convergence threshold
  - `residuals::Vector`, the value of the convergence criteria at each iteration

Note: No warning or error is thrown when the solution is not converged. The user is expected to check convergence manually via the `.isconverged` field.

### Orthogonalization

Orthogonalizing the basis vectors for a `KrylovSubspace` is crucial for numerical stability, and is a core low-level operation for many iterative solvers.

    orthogonalize{T}(v::Vector{T}, K::KrylovSubspace{T}, [p::Int]; [method::Symbol], [normalize::Bool])

- `v`: The vector to orthogonalize
- `K`: The `KrylovSubspace` to orthogonalize against
- `p`: The number of Krylov vectors to orthogonalize against (the default is all available)
- `method`: Orthogonalization method. Currently supported methods are:
  - `:GramSchmidt`. Not recommended, as it is numerically unstable.
  - `:ModifiedGramSchmidt`. (default)
  - `:Householder`.
