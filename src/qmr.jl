export qmr, qmr!

using LinearAlgebra

"""
    qmr(A, b; kwargs...) -> x, [history]

Same as [`qmr!`](@ref), but allocates a solution vector `x` initialized with zeros.
"""
qmr(A, b; kwargs...) = qmr!(zerox(A, b), A, b; initially_zero = true, kwargs...)

"""
    qmr!(x, A, b; kwargs...) -> x, [history]

Solves the problem ``Ax = b`` with the Quasi-Minimal Residual method without lookahead.

# Arguments

- `A`: linear operator;
- `b`: right-hand side.
- `x`: Initial guess, will be updated in-place;

## Keywords

- `initially_zero::Bool`: If `true` assumes that `iszero(x)` so that one
  matrix-vector product can be saved when computing the initial
  residual vector;
- `tol`: relative tolerance;
- `Pl`: left preconditioner;
- `Pr`: right preconditioner;
- `maxiter::Int = size(A, 2)`: Maximum number of iterations
- `log::Bool`: keep track of the residual norm in each iteration;
- `verbose::Bool`: print convergence information during the iterations.

# Return values

**if `log` is `false`**

- `x`: approximate solution.

**if `log` is `true`**

- `x`: approximate solution;
- `history`: convergence history.
"""
function qmr!(x,A, b;
  Pl = I,
  Pr = I,
  tol = sqrt(eps(real(eltype(b)))),
  maxiter = size(A, 2),
  log::Bool = false,
  initially_zero::Bool = false,
  verbose::Bool = false)

  # Startup history file
  history = ConvergenceHistory(partial = !log)
  history[:tol] = tol
  reserve!(history, :resnorm, maxiter)

  # Empty initial videos
  δ = []
  ν = []
  w = []
  y_vec = []
  z_vec = []
  p = []
  q = []
  ϵ = []
  p_vec = []
  β = []
  θ = []

  # Set up initial conditions
  x_last = copy(x)
  r = [b-A*x_last]
  ν_vec = [r[1]]
  y = Pl \ ν_vec[1]
  ρ = [norm(y,2)]
  w_vec = [r[1]]
  z =  Pr' \ w_vec[1]
  ζ = [norm(z,2)]

  γ = [1]
  η = [-1]

  for i in 1:maxiter
    # Check initial state
    if ρ[i] == 0 || ζ[1] == 0
      history.isconverged = false;
      break
    end

    # Generate Lanczos Vectors ν and w
    ν[i] = ν_vec[i]/ρ[i];y /= ρ[i]
    w[i] = w_vec[i]/ζ[i];z /= ζ[i]

    # Another failure case
    δ[i] = dot(z,y)
    if δ[i] == 0
      history.isconverged = false;
      break
    end

    y_vec = Pr \ y
    z_vec = Pl' \ z

    if i == 1
      p[i] = y_vec; q[i] = z_vec
    else
      p[i] = y_vec - ((ζ[i]*δ[i])/(ϵ[i-1]))*p[i-1]
      q[i] = z_vec - ((ρ[i]*δ[i])/(ϵ[i-1]))*q[i-1]
    end

    p_vec = A*p[i]

    # More failure cases
    ϵ[i] = dot(q[i],p_vec)
    β[i] = ϵ[i] / δ[i]
    if ϵ[i] == 0
      history.isconverged = false;
      break;
    end
    if β[i] == 0
      history.isconverged = false;
      break;
    end

    ν_vec[i+1] = p_vec - β[i]*ν[i]
    y = Pl \ ν_vec[i+1]
    ρ[i+1] = norm(y,2)
    w_vec[i+1] = A'*q[i] - β[i]*w[i]
    z = Pr' \ w_vec[i+1]
    ζ[i+1] = norm(z,2)

    # Update QR Factorization
    θ[i] = ρ[i+1] / (γ[i-1] * abs(β[i]))
    γ[i] = 1 / sqrt(1 + θ[i]^2)

    if γ[i] == 0
      history.isconverged = false;
      break;
    end

    η[i] = (-η[i-1]*ρ[i]*γ[i]^2)/(β[i]*γ[i-1]^2)

    if i == 1
      d[i] = η[i]*p[i]; s[i] = η[i]*p_vec
    else
      d[i] = η[i]*p[i] + (θ[i-1]*γ[i])^2*d[i-1]
      s[i] = η[i]*p_vec + (θ[i-1]*γ[i])^2*s[i-1]
    end

    # Update solution and test for convergence
    x[i] = x[i-1] + d[i]
    r[i] = r[i-1] - s[i]
    push!(history, :resnorm, r[i])
    if r[i] <= tol
      history.isconverged = true;
      break;
    end
  end

  log ? (x,history) :  x
end
