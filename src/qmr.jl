import Base: iterate
export qmr, qmr!


# using IterativeSolvers
using LinearAlgebra
using Test
using SparseArrays
using Random

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
function qmr!(x, A, b;
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
  reserve!(history,:resnorm,maxiter)

  # Lambda for right preconditioner
  M1m1x = x-> Pr \x
  M1tm1x = x-> Pr' \x

  # Lamda for left preconditioner
  M2m1x = x-> Pl \x
  M2tm1x = x-> Pl' \x

  # Lambda for multiplying by A
  Ax = x -> A*x
  Atx = x -> A'*x

  r = b - A*x
  bnorm = norm(b)
  res₀ = norm(r)
  res_vec = [res₀]
  vt = r
  y = M1m1x(vt)
  ρ₀ = norm(y)
  wt = r
  z = M2tm1x(wt)
  xi1 = norm(z)
  γ₀ = 1.
  η₀ = -1.

  # Setting up working variables

  # Type of elements in b
  btype = eltype(b)
  ϵ₀ = 0
  v = zeros(btype,length(b))
  w = zeros(btype,length(b))
  p = zeros(btype,length(b))
  q = zeros(btype,length(b))
  β₁ = zero(btype)
  ρ₁ = zero(btype)
  θ₀ = zero(btype)
  θ₁ = zero(btype)
  γ₁ = zero(btype)
  η₁ = zero(btype)
  d = zeros(btype,length(b))
  s = zeros(btype,length(b))
  res₁ = zero(btype)
  last_iter = zero(Int64)

  for iter = 1:maxiter
    ## If ρ₀ == 0 or xi1 == 0, method fails.
    v = vt / ρ₀
    y /= ρ₀
    w = wt / xi1
    z /= xi1
    δ₁ = z' * y   # If δ₁ == 0, method fails.
    yt = M2m1x(y)
    zt = M1tm1x(z)
    if (iter == 1)
      p = yt
      q = zt
    else# 2nd or higher
      p = yt - (xi1*δ₁/ϵ₀) * p
      q = zt - (ρ₀*δ₁/ϵ₀) * q
    end
    pt = Ax(p)
    ϵ₀ = (q' * pt)[1]          # If ϵ₀ == 0, method fails.
    β₁ = ϵ₀ / δ₁   # If β₁ == 0, method fails.
    vt = pt - β₁ * v
    y = M1m1x(vt)
    ρ₁ = norm(y)
    wt = Atx(q) - β₁ * w
    z = M2tm1x(wt)
    xi1 = norm(z)
    θ₁ = ρ₁ / (γ₀ * abs(β₁))
    γ₁ = 1 / sqrt(1 + θ₁^2)   # If γ₁ == 0, method fails.
    η₁ = -η₀ * ρ₀ * γ₁^2 / (β₁ * γ₀^2)

    if (iter == 1)
      d = η₁ * p
      s = η₁ * pt
    else
      d = η₁ * p + (θ₀*γ₁)^2 * d
      s = η₁ * pt + (θ₀ * γ₁)^2 * s
    end
    x += d
    r -= s

    res₁ = norm(r) / bnorm
    push!(res_vec,norm(r))

    # Check for convergance
    if (res₁ < tol)
      # Solver Converged with tolerance
      if log
        history.isconverged = true
      end
      last_iter = iter
      break
    elseif (res₀ <= res₁)
      # Local minimum found
      if log
        history.isconverged = true
      end
      last_iter = iter
      break
    end
    θ₀ = θ₁
    η₀ = η₁
    γ₀ = γ₁
    ρ₀ = ρ₁
    nextiter!(history)
    log && push!(history,:resnorm,norm(r))
  end

  verbose && println()
  log && shrink!(history)

  relres = res₁
  log ? (x, history) :  x
end
