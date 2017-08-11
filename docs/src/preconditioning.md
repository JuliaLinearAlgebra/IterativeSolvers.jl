# [Preconditioning](@id Preconditioning)

Many iterative solvers have the option to provide left and right preconditioners (`Pl` and `Pr` resp.) in order to speed up convergence or prevent stagnation. They transform a problem $Ax = b$ into a better conditioned system $(P_l^{-1}AP_r^{-1})y = P_l^{-1}b$, where $x = P_r^{-1}y$.

These preconditioners should support the operations 

- `A_ldiv_B!(y, P, x)` computes `P \ x` in-place of `y`;
- `A_ldiv_B!(P, x)` computes `P \ x` in-place of `x`;
- and `P \ x`.

If no preconditioners are passed to the solver, the method will default to 

```julia
Pl = Pr = IterativeSolvers.Identity()
```

## Available preconditioners
IterativeSolvers.jl itself does not provide any other preconditioners besides `Identity()`, but recommends the following external packages:

- [ILU.jl](https://github.com/haampie/ILU.jl) for incomplete LU decompositions (using drop tolerance);
- [IncompleteSelectedInversion.jl](https://github.com/ettersi/IncompleteSelectedInversion.jl) for incomplete LDLt decompositions.