# Contributing

Contributions are always welcome, as are feature requests and suggestions. Feel
free to open an [issue](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl/issues) or [pull request](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl/pulls) at any time.

It is important to note that almost every method in the package has documentation,
to know what it does simply use `?<method>` in the terminal.

```julia
julia> using IterativeSolvers

help?> IterativeSolvers.Adivtype
  Adivtype(A, b)

  Determine type of the division of an element of b against an element of A:

  typeof(one(eltype(b))/one(eltype(A)))
```

## Setting workspace up

Julia's internal package manager makes it easy to install and modify packages
from Github. Any package hosted on Github can be installed via `Pkg.clone` by
providing the [repository's URL](https://help.github.com/articles/which-remote-url-should-i-use/), so installing a fork on your system is a simple task.

```julia
Pkg.clone("https://github.com/johndoe/IterativeSolvers.jl")
```

It is to note here if you have the original package installed the fork will
replace it, this is not a problem.

Now find your fork's location.

```julia
Pkg.dir("IterativeSolvers")
```

Once there you will notice you are on the master branch, whenever a package is
imported Julia will use the code in the current branch, this means checking out
other git branches will let you use/test whatever there is.

## Adding or modifying iterative methods

Each iterative method method must log information using the inner `ConvergenceHistory` type.

There are two types of `ConvergenceHistory`: plain and restarted. The only real
difference between the two is how they are plotted and how the number of restarts
is calculated, everything else is the same.

Before logging information space must always be reserved.

```julia
log = ConvergenceHistory()
log[:tol] = tol
reserve!(log,:betas, maxiter) # Vector of length maxiter
reserve!(log,:conv, maxiter, T=BitArray) # Vector of length maxiter
reserve!(log,:ritz, maxiter, k) # Matrix of size (maxiter, k)
```

To store information at each iteration use `push!`.

```julia
push!(log, :conv, conv)
push!(log, :ritz, F[:S][1:k])
push!(log, :betas, L.β)
```

To advance the log index to the next iteration use `nextiter!`.

```julia
nextiter!(log)
```

A more detailed explanation of all the functions is in both the public and internal
documentation of `ConvergenceHistory`.

The most rich example of the usage of `ConvergenceHistory` is in `svdl`.
