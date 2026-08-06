# [LALQMR](@id LALQMR)

Look-ahead Lanczos Quasi-minimal Residual (LALQMR) is the look-ahead variant of [QMR](@ref) for solving $Ax = b$ approximately for $x$ where $A$ is a linear operator and $b$ the right-hand side vector. $A$ may be non-symmetric [^Freund1990]. The Krylov subspace is generated via the [Look-ahead Lanczos process](@ref LAL)

## Usage

```@docs
lalqmr
lalqmr!
```

## Implementation details
This implementation of LALQMR follows [^Freund1994] and [^Freund1993], where we generate the Krylov subspace via Lanczos bi-orthogonalization based on the matrix `A` and its transpose. In the regular Lanczos process, we may encounter singularities during the construction of the Krylov basis. The look-ahead process avoids this by building blocks instead, avoiding the singularities. Therefore, the LALQMR technique is well-suited towards problematic linear systems. Typically, most systems will not have blocks of size larger than 4 made, and most iterations are "regular", or blocks of size 1 [^Freund1994].

For more detail on the implementation see the original paper [^Freund1994]. This implementation follows the paper, with the exception that if a block of maximum size is reached during the Lanczos process, the block is closed instead of restarted.

For a solution vector of size `n`, the memory allocated will be `(12 + max_memory * 5) * n`.

!!! tip
    LALQMR can be used as an [iterator](@ref Iterators) via `lalqmr_iterable!`. This makes it possible to access the next, current, and previous Krylov basis vectors during the iteration.

## References
[^Freund1993]:
Freund, R. W., Gutknecht, M. H., & Nachtigal, N. M. (1993). An Implementation of the Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices. SIAM Journal on Scientific Computing, 14(1), 137–158. https://doi.org/10.1137/0914009

[^Freund1994]:
Freund, R. W., & Nachtigal, N. M. (1994). An Implementation of the QMR Method Based on Coupled Two-Term Recurrences. SIAM Journal on Scientific Computing, 15(2), 313–337. https://doi.org/10.1137/0915022
[^Freund1990]:

Freund, W. R., & Nachtigal, N. M. (1990). QMR : for a Quasi-Minimal Residual Linear Method Systems. (December).
