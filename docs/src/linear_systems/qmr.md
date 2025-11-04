# [QMR](@id QMR)

QMR is a short-recurrence version of GMRES for solving $Ax = b$ approximately for $x$ where $A$ is a linear operator and $b$ the right-hand side vector. $A$ may be non-symmetric.

## Usage

```@docs
qmr
qmr!
```

## Implementation details
QMR exploits the tridiagonal structure of the Hessenberg matrix. Although QMR is similar to GMRES, where instead of using the Arnoldi process, a pair of biorthogonal vector spaces $V$ and $W$ is constructed via the Lanczos process. It requires that the adjoint of $A$ `adjoint(A)` be available [^Saad2003], [^Freund1990].

QMR enables the computation of $V$ and $W$ via a three-term recurrence. A three-term recurrence for the projection onto the solution vector can also be constructed from these values, using the portion of the last column of the Hessenberg matrix. Therefore we pre-allocate only eight vectors.

!!! tip
    QMR can be used as an [iterator](@ref Iterators) via `qmr_iterable!`. This makes it possible to access the next, current, and previous Krylov basis vectors during the iteration.
