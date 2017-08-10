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
vectors to vectors. Typically `A` is a `Matrix` or a `SparseMatrixCSC`, but since
`A` is not explicitly typed, any linear operator that supports matrix operations
can be used as well. This makes it possible to apply solvers *matrix-free*. In 
IterativeSolvers.jl we strongly recommend [LinearMaps.jl](https://github.com/Jutho/LinearMaps.jl) 
for non-matrix types of `A`.

For matrix-free types of `A` the following interface is expected to be defined:

- `A*v` computes the matrix-vector product on a `v::AbstractVector`;
- `A_mul_B!(y, A, v)` computes the matrix-vector product on a `v::AbstractVector` in-place;
- `eltype(A)` returns the element type implicit in the equivalent matrix representation of `A`;
- `size(A, d)` returns the nominal dimensions along the `d`th axis in the equivalent matrix representation of `A`.

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

`A` is a linear operator as described above.

`b` is the vector to be solved.

`x` is a vector for the initial guess. In the case of a mutating call this
parameter will be overwritten.

Output will be the solution to the system.


### Additional arguments

Keyword names will vary depending on the method, however some of them will always have the same spelling:

- `tol`: (relative) stopping tolerance of the method;
- `verbose`: print information during the iterations;
- `maxiter`: maximum number of allowed iterations;
- `Pl` and `Pr`: left and right preconditioner. See [Preconditioning](@ref);
- `log::Bool = false`: output an extra element of type `ConvergenceHistory` containing the convergence history.

### `log` keyword

Most solvers contain the `log` keyword. This is to be used when obtaining
more information is required, to use it place the set `log` to `true`.

```julia
x, ch = cg(Master, rand(10, 10), rand(10) log=true)
svd, L, ch = svdl(Master, rand(100, 100), log=true)
```

The function will now return one more parameter of type `ConvergenceHistory`.

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

`sep::Symbol = :white`: color of the line separator in restarted methods.
