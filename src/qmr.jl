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
function qmr!(x,A, b;
  Pl = I,
  Pr = I,
  tol = sqrt(eps(real(eltype(b)))),
  maxiter = size(A, 2),
  log::Bool = false,
  initially_zero::Bool = false,
  verbose::Bool = false)
  # Implemented using qmr method given in
  #  [1] https://link.springer.com/content/pdf/10.1007%2FBF01385726.pdf
  # Algorithm 3.1

  # Test if singular then break
  # Test if Hermitian -> use conjugate gradient

  # Startup history file
  history = ConvergenceHistory(partial = !log)
  history[:tol] = tol
  reserve!(history, :resnorm, maxiter)

  # Empty initial vectors
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
  d = []
  s = []

  # [1] 3.1.0

  # Choose x₀ ∈ Cᴺ and set r₀ = b - Ax₀, ρ₀ = |r₀|, v₁ = r₀/ρ₀
  # Choose w₁ ∈ Cᴺ with |w₁|= 1
  # Set up initial conditions
  x_vec = [copy(x)]
  r = [b-A*x_vec[1]]
  ν_vec = [r[1]]
  y = Pl \ ν_vec[1] #? Why solving within solver
  ρ = [norm(y,2)]
  w_vec = [r[1]]
  z =  Pr' \ w_vec[1]#? Why solving within solver
  ζ = [norm(z,2)]

  γ = [1.]
  η = [-1.]

  for i in 1:maxiter
    # [1] 3.1.1 Perform the nth iteration of the look-ahead Lanczos Algorithm 2.1;
    # This yields matrices Vⁿ, Vⁿ⁺¹, Hⁿₑ which satisfy (3.5);
    # Check initial state
    if ρ[i] == 0 || ζ[1] == 0
      history.isconverged = false
      break
    end

    # Generate Lanczos Vectors ν and w
    push!(ν, ν_vec[i]/ρ[i])
    y /= ρ[i]
    push!(w,w_vec[i]/ζ[i])
    z /= ζ[i]

    # Another failure case
    push!(δ,dot(z,y))
    if δ[i] == 0
      history.isconverged = false
      break
    end

    y_vec = Pr \ y
    z_vec = Pl' \ z

    if i == 1
      push!(p,y_vec)
      push!(q,z_vec)
    else
      push!(p,y_vec - ((ζ[i]*δ[i])/(ϵ[i-1]))*p[i-1])
      push!(q,z_vec - ((ρ[i]*δ[i])/(ϵ[i-1]))*q[i-1])
    end

    p_vec = A*p[i]

    # More failure cases
    push!(ϵ,dot(q[i],p_vec))
    push!(β,ϵ[i] / δ[i])
    if ϵ[i] == 0
      history.isconverged = false
      break
    end
    if β[i] == 0
      history.isconverged = false
      break
    end

    push!(ν_vec,p_vec - β[i]*ν[i])
    y = Pl \ ν_vec[i+1]
    push!(ρ,norm(y,2))
    push!(w_vec,A'*q[i] - β[i]*w[i])
    z = Pr' \ w_vec[i+1]
    push!(ζ,norm(z,2))

    # Update QR Factorization
    push!(θ, ρ[i+1] / (γ[i] * abs(β[i])))
    push!(γ,1 / sqrt(1 + θ[i]^2))

    if γ[i+1] == 0
      history.isconverged = false
      break
    end

    push!(η,(-η[i]*ρ[i]*γ[i]^2)/(β[i]*γ[i]^2))

    if i == 1
      push!(d,η[i]*p[i])
      push!(s,η[i]*p_vec)
    else
      push!(d,η[i]*p[i] + (θ[i-1]*γ[i])^2*d[i-1])
      push!(s,η[i]*p_vec + (θ[i-1]*γ[i])^2*s[i-1])
    end

    # Update solution and test for convergence
    push!(x_vec,x_vec[i] + d[i])
    push!(r,r[i] - s[i])
    push!(history, :resnorm, norm(r[i+1]))
    if norm(r[i+1]) <= tol
      history.isconverged = true
      break
    end
  end

  log ? (x_vec[end],history) :  x_vec[end]
end


# Testing
n = 4
T = ComplexF64
A = spdiagm(-1 => fill(-1.0,n-1), 1 => fill(4.0,n-1))
b = sum(A,dims=2)
M1 = spdiagm(-1 => fill(-1/2,n-1), 0 => ones(n))
M2 = spdiagm(0 => fill(4.0 ,n), 1 => fill(-1.0,n-1))
# x = ones(T,100)
x0 = rand(T,n)
x = qmr!(x0,A,b)

inv(Array(A))*b
