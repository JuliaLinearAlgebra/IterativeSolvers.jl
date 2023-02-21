# Iterative Solvers

Fork 自[IterativeSolvers.jl](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl)。存在以下改动：
1. Verbose 动态输出增加相对残差列，使收敛速度更直观；
2. gmres 求解器的基向量空间允许 AbstractArray 以方便移植。

[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl/blob/master/LICENSE)
[![Build Status](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl/workflows/CI/badge.svg)](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl/actions)
[![Codecov](http://codecov.io/github/JuliaLinearAlgebra/IterativeSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLinearAlgebra/IterativeSolvers.jl?branch=master)
[![Coveralls](https://coveralls.io/repos/JuliaLinearAlgebra/IterativeSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaLinearAlgebra/IterativeSolvers.jl?branch=master)

`IterativeSolvers` is a [Julia](http://julialang.org) package that provides iterative algorithms for solving linear systems, eigensystems, and singular value problems.

## Resources

- [Manual](https://IterativeSolvers.julialinearalgebra.org/dev/)
- [Contributing](https://julialinearalgebra.github.io/IterativeSolvers.jl/dev/about/CONTRIBUTING/)

## Installing

To install the package, open the package manager in the REPL via `]` and run

```julia
pkg> add IterativeSolvers
```
