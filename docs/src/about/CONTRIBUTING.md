# Contributing

Contributions are always welcome, as are feature requests and suggestions. Feel
free to open [issues](https://help.github.com/articles/creating-an-issue/) and [pull requests](https://help.github.com/articles/creating-a-pull-request/) at any time. If you aren't familiar with git or Github please start [now](https://help.github.com/articles/good-resources-for-learning-git-and-github/).

## Setting workspace up

Julia's internal package manager makes it easy to install and modify packages
from Github. Any package hosted on Github can be installed via `Pkg.clone` by
providing the [repository's URL](https://help.github.com/articles/which-remote-url-should-i-use/), so installing a fork on your system is a simple task.

```julia
Pkg.clone("https://github.com/johndoe/IterativeSolvers.jl")
```

It is to note here if you have the original package installed the fork will
replace it, this is not a problem.

Now find your fork's location

```julia
Pkg.dir("IterativeSolvers")
```

Once there you will notice you are on the master branch, whenever a package is
imported Julia will use the code in the current branch, this means checking out
other git branches will let you use/test whatever there is.

## Adding benchmarks

The [Benchmarks](@ref) tab of the documentation is built automatically with Travis.
Any benchmark added will be displayed automatically after a successful pull request.

The benchmark suite gets built doing a cross product between the available matrices
and available methods, if there are `n` methods and `m` linear operators then `n*m` will be
the upper limit of benchmarks to be made. Some methods are not compatible with certain
matrices, to avoid generating unnecessary benchmarks each method and matrix has
traits, linear operator traits are inspired from [MatrixDepot.jl](http://matrixdepotjl.readthedocs.io/en/latest/properties.html).

**Method traits**

* accessible : Method accesses the linear operator's fields.
* inverse    : `A`'s Inverse must exist.
* symmetric  : `A`'s must be symmetric.
* pos-def    : `A`'s must be definite.

**Linear Operator traits**

* accessible : Is accessible.
* inverse    : `A` is exist.
* symmetric  : `A` is symmetric.
* pos-def    : `A` is definite.
* eigen      : Part of the eigensystem of the matrix is explicitly known.
* graph      : An adjacency matrix of a graph.
* ill-cond   : The matrix is ill-conditioned for some parameter values.
* random     : The matrix has random entries.
* regprob    : The output is a test problem for Regularization Methods.
* sparse     : The matrix is sparse.

A benchmark between a method and a linear operator will be made if and only if
the traits of the method is subset of the traits of the linear operator.

Benchmarks are stored in [Benchmarks.jl](https://github.com/JuliaLang/IterativeSolvers.jl/tree/master/benchmark/Benchmarks.jl).
To add a method use `addEqMethod`.

```julia
addEqMethod(methods, "jacobi", jacobi, ["inverse","accessible"])
addEqMethod(methods, "gauss_seidel", gauss_seidel, ["inverse","accessible"])
addEqMethod(methods, "sor", sor, ["inverse","accessible"])
addEqMethod(methods, "ssor", ssor, ["inverse","accessible", "symmetric"])
addEqMethod(methods, "cg", cg, ["inverse", "symmetric", "pos-def"])
addEqMethod(methods, "gmres", gmres, ["inverse"])
addEqMethod(methods, "lsqr", lsqr, ["inverse"])
addEqMethod(methods, "chebyshev", chebyshev, ["inverse", "accessible"])
```

Here `methods` is a dictionary, the second argument is the name to be displayed in
the benchmarks, the third argument is the function and the fourth is the traits.
Every function has a predetermined call in `buildCall` function.

To add an equation use `addEquation`.

```julia
#Sparse matrix equations
addEquation(
    equations, "Poisson", ["Sparse", "Poisson"],
    ["sparse","inverse", "symmetric", "pos-def", "eigen", "accessible"],
    :(matrixdepot("poisson",4))
    )

#Function matrix equations
addEquation(
    equations, "SOLtest", ["Function", "SOLtest"],
    ["function","inverse"],
    :(buildSol(10)),
    10
)
```

Here `equations` is a dictionary, the second argument is the name to be displayed in
the benchmarks, the third argument is the path inside the `BenchmarkGroup` type
the fourth argument is the traits, the fifth is the matrix generator and
the sixth is the size of the matrix. The size of the matrix has to be passed when
it is impossible to deduce the dimension from the generator, in this case buildSol
generates a function and not a matrix.

To add a custom benchmark use directly the `suite` variable which is the `BenchmarkGroup`
of the package, to know more of this type check [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).
