# [Look-Ahead Lanczos Bi-Orthogonalization](@id LAL)

The Look-ahead Lanczos Bi-orthogonalization process is a generalization of the Lanczos
bi-orthoganilization (as used in, for instance, [QMR](@ref)). The look-ahead process detects
(near-)singularities during the construction of the Lanczos iterates, known as break-down,
and skips over them. This is particularly advantageous for linear systems with large
null-spaces or eigenvalues close to 0.

We provide an iterator interface for the Lanczos decomposition, following [^Freund1994],
where the look-ahead Lanczos process is implemented as a two-term coupled recurrence.

## Usage

```@docs
LookAheadLanczosDecomp
```

## References

[^Freund1993]: Freund, R. W., Gutknecht, M. H., & Nachtigal, N. M. (1993). An Implementation
    of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices. SIAM Journal on
    Scientific Computing, 14(1), 137–158. https://doi.org/10.1137/0914009

[^Freund1994]: Freund, R. W., & Nachtigal, N. M. (1994). An Implementation of the QMR Method
    Based on Coupled Two-Term Recurrences. SIAM Journal on Scientific Computing, 15(2),
    313–337. https://doi.org/10.1137/0915022