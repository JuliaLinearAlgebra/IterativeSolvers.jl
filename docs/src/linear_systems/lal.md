# [Look-Ahead Lanczos Bi-Orthogonalization](@id LAL)

The Look-ahead Lanczos Bi-orthogonalization process is a generalization of the Lanczos
bi-orthoganilization (as used in, for instance, [QMR](@ref)) [^Freund1993], [^Freund1994]. The look-ahead process detects
(near-)singularities during the construction of the Lanczos iterates, known as break-down,
and skips over them. This is particularly advantageous for linear systems with large
null-spaces or eigenvalues close to 0.

We provide an iterator interface for the Lanczos decomposition (Freund, 1994),
where the look-ahead Lanczos process is implemented as a two-term coupled recurrence.

## Usage

```@docs
IterativeSolvers.LookAheadLanczosDecomp
IterativeSolvers.LookAheadLanczosDecompLog
```
