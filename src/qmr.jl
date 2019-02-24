export qmr, qmr!

using LinearAlgebra

"""
    qmr(A, b; kwargs...) -> x, [history]

Same as [`qmr!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
qmr(A, b; kwargs...) = qmr!(zerox(A, b), A, b; initially_zero = true, kwargs...)

"""
    qmr!(x, A, b; kwargs...) -> x, [history]

Solves the problem ``Ax = b`` with the Quasi-Minimal Residual method.

# Arguments

<<<<<<< HEAD
- `A`: linear operator;
- `b`: right-hand side.
- `x`: Initial guess, will be updated in-place;
=======
- `x`: Initial guess, will be updated in-place;
- `A`: linear operator;
- `b`: right-hand side.
>>>>>>> 8c9ae8a8aee84bfa4974c21cfb134c11ab265870

## Keywords

- `initially_zero::Bool`: If `true` assumes that `iszero(x)` so that one
  matrix-vector product can be saved when computing the initial
  residual vector;
- `tol`: relative tolerance;
- `Pl`: left preconditioner;
- `Pr`: right preconditioner;
<<<<<<< HEAD
- `maxiter::Int = size(A, 2)`: Maximum number of iterations
=======
>>>>>>> 8c9ae8a8aee84bfa4974c21cfb134c11ab265870
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iterations.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;
- `history`: convergence history.
"""
<<<<<<< HEAD

#########################
# Method Implementation #
#########################

# Kiran Shila and Eric Valentino
# Dept of Electrical Engineering
# University of South Florida
#-----------------------------
function qmr!(x,A, b;
  Pl = Identity(),
  Pr = Identity(),
  tol = sqrt(eps(real(eltype(b)))),
  maxiter = size(A, 2),
=======
function qmr!(x, A, b;
  Pl = Identity(),
  Pr = Identity(),
  tol = sqrt(eps(real(eltype(b)))),
>>>>>>> 8c9ae8a8aee84bfa4974c21cfb134c11ab265870
  log::Bool = false,
  initially_zero::Bool = false,
  verbose::Bool = false)


end
