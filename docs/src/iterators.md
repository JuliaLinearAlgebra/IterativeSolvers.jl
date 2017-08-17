# [Iterative solvers as iterators](@id Iterators)

In advanced use cases you might want to access the internal data structures of the solver, inject code to be run after each iteration, have total control over allocations or reduce overhead in initialization. The iterator approach of IterativeSolvers.jl makes this possible.

!!! note
    At this point [BiCGStab(l)](@ref BiCGStabl), [CG](@ref CG), [Chebyshev](@ref Chebyshev), [GMRES](@ref GMRES), [MINRES](@ref MINRES) and the [stationary methods](@ref Stationary) are implemented as iterators. However, the package does not yet export the iterators and helper methods themselves.

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

```julia
A = sprand(10_000, 10_000, 10 / 10_000) + 20I
b1 = rand(10_000)
b2 = rand(10_000)
x = rand(10_000)

my_iterable = IterativeSolvers.jacobi_iterable(x, A, b1, maxiter = 4)

for item in my_iterable 
    println("Iteration for rhs 1")
end

@show norm(b1 - A * x) / norm(b1)

# Copy the next right-hand side into the iterable
copy!(my_iterable.b, b2)

for item in my_iterable
    println("Iteration for rhs 2")
end

@show norm(b2 - A * x) / norm(b2)
```

This would output:

```
Iteration for rhs 1
Iteration for rhs 1
Iteration for rhs 1
Iteration for rhs 1
norm(b1 - A * x) / norm(b1) = 0.08388528015119746
Iteration for rhs 2
Iteration for rhs 2
Iteration for rhs 2
Iteration for rhs 2
norm(b2 - A * x) / norm(b2) = 0.0003681972775644809
```

## Other use cases
Other use cases include: 
- computing the (harmonic) Ritz values from the Hessenberg matrix in GMRES;
- comparing the approximate residual of methods such as GMRES and BiCGStab(l) with the true residual during the iterations;
- updating a preconditioner in flexible methods.

