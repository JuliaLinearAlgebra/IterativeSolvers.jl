# Getting started

## Installation

The package can be installed via Julia's package manager.

```julia
julia> Pkg.add("IterativeSolvers")
```

## Interface

Virtually all solvers have the common function declarations:

```julia
solver(A, args...; kwargs...)
solver!(x, A, args...; kwargs...)
```

where `A` is a [linear operator](@ref matrixfree) and `x` an initial guess. The second declaration also updates `x` in-place.

### [Explicit matrices and the matrix-free approach](@id matrixfree)
Rather than constructing an explicit matrix `A` of the type `Matrix` or `SparseMatrixCSC`, it is also possible to pass a general linear operator that performs matrix operations implicitly. This is called the **matrix-free** approach.

For matrix-free types of `A` the following interface is expected to be defined:

- `A*v` computes the matrix-vector product on a `v::AbstractVector`;
- `A_mul_B!(y, A, v)` computes the matrix-vector product on a `v::AbstractVector` in-place;
- `eltype(A)` returns the element type implicit in the equivalent matrix representation of `A`;
- `size(A, d)` returns the nominal dimensions along the `d`th axis in the equivalent matrix representation of `A`.

!!! tip "Matrix-free with LinearMaps.jl"
    We strongly recommend [LinearMaps.jl](https://github.com/Jutho/LinearMaps.jl) for matrix-free linear operators, as it implements the above methods already for you; you just have to write the action of the linear map.


### Additional arguments

Keyword names will vary depending on the method, however some of them will always have the same spelling:

- `tol`: (relative) stopping tolerance of the method;
- `verbose`: print information during the iterations;
- `maxiter`: maximum number of allowed iterations;
- `Pl` and `Pr`: left and right preconditioner. See [Preconditioning](@ref Preconditioning);
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

```@docs
ConvergenceHistory
```
