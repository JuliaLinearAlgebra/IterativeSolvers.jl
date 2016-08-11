# Manual

## Installation

The package can be installed with a simple instruction.

```julia
julia> Pkg.add("IterativeSolvers")
```

After installing the package, if you wish to use the latest features of the
package you must switch to the master branch with `Pkg.clone`.

```julia
julia> Pkg.checkout("IterativeSolvers")
```

## Interface

All linear-algebraic routines will take as input a linear operator `A` that maps
vectors to vectors. `A` is not explicitly typed, but must either be a
[`KrylovSubspace`](@ref) or support multiplication `*` or function composition (apply)
that behave as necessary to produce the correct mapping on the vector space.

A custom type for `A` may be specified. The following interface is expected to
be defined on `A`:

* `A*v` is defined and computes the matrix-vector product on a `v::Vector`.

* `eltype(A)` is defined and returns the element type implicit in the equivalent
matrix representation of `A`.

* `size(A, d)` is defined and returns the nominal dimensions along the dth axis
in the equivalent matrix representation of `A`.

### Solvers

All linear solvers have a common function declaration (with a few exceptions).

```
solver(A, b::Vector; kwargs...)
solver!(x, A, b::Vector; kwargs...)
```

In the case of eigenproblems or singular value decompositions:

```
eigsolver(A; kwargs...)
eigsolver!(x, A; kwargs...)
```

* `A` is a linear operator as described above.

* `b` is the vector to be solved.

* `x` is a vector for the initial guess. In the case of a mutating call this
parameter will be overwritten. (Mutating functions end with `!`)

* output will be the solution the system.


### Additional arguments

Keyword names will vary depending on the method, however some
of them will always have the same spelling:

* `tol`: stopping tolerance of the method. When a method accepts more than one
tolerance they are enumerated  with a letter prefix, e.g `atol`, `btol`, `ctol`,
etc.

* `verbose`: print information about the running method.

* `maxiter`: maximum number of allowed iterations.

* `Pl`: left preconditioner. (When applicable)

* `Pr`: right preconditioner. (When applicable)

* `log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

* `plot`: plot information relevant to the method. (Only for `Master` version)

### `log` keyword

All solvers contain the `log` keyword. This is to be used when obtaining
more information is required, to use it place the set `log` to `true`.

```julia
x, ch = cg(Master, rand(10, 10), rand(10) log=true)
svd, L, ch = svdl(Master, rand(100, 100), plot=true, log=true)
```

The function will now return one more parameter of type [`ConvergenceHistory`](@ref).

*`Note:`*  Keyword argument `plot` is only available when `log` is set.

## ConvergenceHistory

A [`ConvergenceHistory`](@ref) instance stores information of a solver.

Number of iterations.

```julia
ch.iters
```

Convergence status.

```julia
ch.isconverged
```

Stopping tolerances. (A `Symbol` key is needed to access)

```julia
ch[:tol]
```

Maximum number of iterations per restart. (Only on restarted methods)

```julia
nrests(ch)
```

Number of matrix-vectors and matrix-transposed-vector products.

```julia
nprods(ch)
```

Data stored on each iteration, accessed information can be either a vector
or matrix. This data can be a lot of things, most commonly residual.
(A `Symbol` key is needed to access)

```julia
ch[:resnorm] #Vector or Matrix
ch[:resnorm, x] #Vector or Matrix element
ch[:resnorm, x, y] #Matrix element
```

The available keys of each method is described in the [Public Documentation](@ref).

### Plotting

`ConvergeHistory` provides a recipe to use with the package [Plots.jl](https://github.com/tbreloff/Plots.jl), this makes it really easy to
plot on different plot backends. There are two recipes provided:

One for the whole `ConvergenceHistory`.

```julia
plot(ch)
```

The other one to plot data binded to a key.

```julia
_, ch = gmres(rand(10,10), rand(10), maxiter = 100, log=true)
plot(ch, :resnorm, sep = :blue)
```

*Plot additional keywords*

* `sep::Symbol = :white`: color of the line separator in restarted methods.

## KrylovSubspace

When [`KrylovSubspace`](@ref) is supported by the method, `A` can be replaced
by an instance `K` of it. To check if a certain function supports it use  
`methods` or `?` to find out, if there is a `K` in the argument list then it does.

```julia
julia> ?cg!

    cg!(x, A, b)
    cg!(x, K, b)

    ...
```

`KrylovSubspace` is only allowed on functions that accept `x` as an argument,
that's why `cg` doesn't allow it. This is also true when `x` is a keyword as
is the case of `powm`.

```julia
julia> ?cg!

    powm(A)
    powm(K)

    ...
```

The `KrylovSubspace` type collects information on the Krylov subspace generated
over the course of an iterative Krylov solver.

Recall that the Krylov subspace of order r given a starting vector `b` and a
linear operator `A` is spanned by the vectors `[b, A*b, A^2*b,... A^(r-1)*b]`.
Many modern iterative solvers work on Krylov spaces which expand until they span
enough of the range of `A` for the solution vector to converge. Most practical
algorithms, however, will truncate the order of the Krylov subspace to save
memory and control the accumulation of roundoff errors. Additionally, they do
not work directly with the raw iterates `A^n*b`, but will orthogonalize
subsequent vectors as they are computed so as to improve numerical stability.
`KrylovSubspace`s provide a natural framework to express operations such as a
(possibly non-orthonormalized) basis for the Krylov subspace, retrieving the
next vector in the subspace, and orthogonalizing an arbitrary vector against
(part or all of) the subspace.

The implementation of `KrylovSubspace` in this package differs from standard
textbook presentations of iterative solvers. First, the `KrylovSubspace` type
shows clearly the relationship between the linear operator A and the sequence of
basis vectors for the Krylov subspace that is generated by each new iteration.
Second, the grouping together of basis vectors also makes the orthogonalization
steps in each iteration routine clearer. Third, unlike in other languages, the
powerful type system of Julia allows for a simplified notation for iterative
algorithms without compromising on performance, and enhances code reuse across
multiple solvers.

### Constructors

A [`KrylovSubspace`](@ref) can be initialized by three constructors, depending on the type
of the linear operator:

```
KrylovSubspace{T}(A::AbstractMatrix{T}, [order::Int, v::Vector{Vector{T}}])

KrylovSubspace{T}(A::KrylovSubspace{T}, [order::Int, v::Vector{Vector{T}}])

KrylovSubspace(A, n::Int, [order::Int, v::Vector{Vector{T}}])
```

* `A`: The linear operator associated with the `KrylovSubspace`.

* `order`: the order of the `KrylovSubspace`, i.e. the maximal number of Krylov
vectors to remember.

* `n`: the dimensionality of the underlying vector space `T^n`.

* `v`: the iterable collection of Krylov vectors (of maximal length order).

The dimensionality of the underlying vector space is automatically inferred where
possible, i.e. when the linear operator is an `AbstractMatrix` or `KrylovSubspace`.

*`Note:`* second constructor destroys the old `KrylovSubspace`.

### Orthogonalization

Orthogonalizing the basis vectors for a `KrylovSubspace` is crucial for numerical
stability, and is a core low-level operation for many iterative solvers.

```
orthogonalize{T}(v::Vector{T}, K::KrylovSubspace{T}, [p::Int]; [method::Symbol], [normalize::Bool])
```

* `v`: vector to orthogonalize.

* `K`: `KrylovSubspace` to orthogonalize against.

* `p`: number of Krylov vectors to orthogonalize against. (Default is all available)

* `method`: Orthogonalization method. Currently supported methods are:
`:GramSchmidt`, `:ModifiedGramSchmidt` (default) and `:Householder`.

## Define function as matrix

Suppose you have a function implementing multiplication of `b` by `A`.
This function can have any syntax, but for the purposes of illustration let's
suppose it's defined as:

```julia
mulbyA!(output, b, Adata)
```

Where `Adata` might be some parameters that your function needs.

You can represent it as a linear operator using:

```julia
A = FuncMat(m, n, typ=T, mul = (output, b) -> mulbyA!(output, b, Adata))
```

where `T` is the "element type" of `Adata`.
Note that there are a couple of requirements:

  - `mulbyA!` stores the result in the pre-allocated output.
  - `mulbyA!` should also return output as its sole return value.

If the algorithm also needs multiplication by A', use keyword `cmul`.

```julia
A = FuncMat(
    m, n, typ=T,
    mul = (output, b) -> mulbyA!(output, b, Adata),
    cmul = (output, b) -> mulbyActrans!(output, b, Adata)
    )
```
