# [Preconditioning](@id Preconditioning)

Many iterative solvers have the option to provide left and right preconditioners (`Pl` and `Pr` resp.) in order to speed up convergence or prevent stagnation. They transform a problem $Ax = b$ into a better conditioned system $(P_l^{-1}AP_r^{-1})y = P_l^{-1}b$, where $x = P_r^{-1}y$.

These preconditioners should support the operations

- `ldiv!(y, P, x)` computes `P \ x` in-place of `y`;
- `ldiv!(P, x)` computes `P \ x` in-place of `x`;
- and `P \ x`.

If no preconditioners are passed to the solver, the method will default to

```julia
Pl = Pr = IterativeSolvers.Identity()
```

## Available preconditioners
IterativeSolvers.jl itself does not provide any other preconditioners besides `Identity()`, but recommends the following external packages:

- [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl) for incomplete LU decompositions (using drop tolerance).
- [IncompleteSelectedInversion.jl](https://github.com/ettersi/IncompleteSelectedInversion.jl) for incomplete LDLt decompositions.
- [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) for some algebraic multigrid (AMG) preconditioners.
- [Preconditioners.jl](https://github.com/mohamed82008/Preconditioners.jl) which wraps a bunch of preconditioners from other packages. If you are a beginner or want to try different ones quickly, this is good starting place.
