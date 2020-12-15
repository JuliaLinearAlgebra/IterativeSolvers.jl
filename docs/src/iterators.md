# [Iterative solvers as iterators](@id Iterators)

In advanced use cases you might want to access the internal data structures of the solver, inject code to be run after each iteration, have total control over allocations or reduce overhead in initialization. The iterator approach of IterativeSolvers.jl makes this possible.

!!! note
    At this point [BiCGStab(l)](@ref BiCGStabl), [CG](@ref CG), [Chebyshev](@ref Chebyshev), [GMRES](@ref GMRES), [MINRES](@ref MINRES), [QMR](@ref QMR) and the [stationary methods](@ref Stationary) are implemented as iterators. However, the package does not yet export the iterators and helper methods themselves.

## How iterators are implemented
The solvers listed above are basically a thin wrapper around an iterator. Among other things, they initialize the iterable, loop through the iterator and return the result:

```julia
function my_solver!(x, A, b)
    iterable = MySolverIterable(x, A, b)
    for item in iterable end
    return iterable.x
end
```

Rather than calling `my_solver!(x, A, b)`, you could also initialize the iterable yourself and perform the `for` loop.

## Example: avoiding unnecessary initialization
The Jacobi method for `SparseMatrixCSC` has some overhead in intialization; not only do we need to allocate a temporary vector, we also have to search for indices of the diagonal (and check their values are nonzero). The current implementation initializes the iterable as:

```julia
jacobi_iterable(x, A::SparseMatrixCSC, b; maxiter::Int = 10) =
    JacobiIterable{eltype(x), typeof(x)}(OffDiagonal(A, DiagonalIndices(A)), x, similar(x), b, maxiter)
```

Now if you apply Jacobi iteration multiple times with the same matrix for just a few iterations, it makes sense to initialize the iterable only once and reuse it afterwards:

```jldoctest
julia> using LinearAlgebra, SparseArrays, IterativeSolvers

julia> A = spdiagm(-1 => -ones(3), 0 => 2*ones(4), 1 => -ones(3));

julia> b1 = [1.0, 2, 3, 4];

julia> b2 = [-1.0, 1, -1, 1];

julia> x = [0.0, -1, 1, 0];

julia> my_iterable = IterativeSolvers.jacobi_iterable(x, A, b1, maxiter = 2);

julia> norm(b1 - A * x) / norm(b1)
1.2909944487358056

julia> for item in my_iterable
           println("Iteration for rhs 1")
       end
Iteration for rhs 1
Iteration for rhs 1

julia> norm(b1 - A * x) / norm(b1)
0.8228507357554791

julia> # Copy the next right-hand side into the iterable
       copyto!(my_iterable.b, b2);

julia> norm(b2 - A * x) / norm(b2)
2.6368778887161235

julia> for item in my_iterable
           println("Iteration for rhs 2")
       end
Iteration for rhs 2
Iteration for rhs 2

julia> norm(b2 - A * x) / norm(b2)
1.610815496107484
```


## Other use cases
Other use cases include:
- computing the (harmonic) Ritz values from the Hessenberg matrix in GMRES;
- comparing the approximate residual of methods such as GMRES and BiCGStab(l) with the true residual during the iterations;
- updating a preconditioner in flexible methods.
