using UnicodePlots
import Base: \, eltype, empty!, eps, length, ndims, push!, real, size, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export A_mul_B

abstract Residuals
type ResNone <: Residuals end
type ResArray <: Residuals
  top::Int
  residuals::Array{Float64,1}
end
type RestResArray <: Residuals
  top::Array{Int,1}
  residuals::Array{Float64,2}
  restart::Int
end
Resnorms(maxiter) = ResArray(1,zeros(maxiter))
Resnorms(maxiter,restart) = RestResArray([1,1],zeros(maxiter,restart),restart)

ResArray(dims,maxiter) = Residuals(dims,ones(length(dims)),zeros(dims...))
add!(rn::ResNone) = rn
function add!(ra::ResArray, res::Real)
  ra.residuals[ra.top] = res
  ra.top+=1
  ra
end
function add!(ra::RestResArray, res::Real)
  ra.residuals[ra.top...] = res
  ra.top[2]+=1
  if ra.top[2] > ra.restart
    ra.top[2] = 1
    ra.top[1] += 1
  end
  ra
end
extract(ra::ResArray) = ra.residuals[1:ra.top-1]
function extract(ra::RestResArray)
  i = ra.top[1]
  ra.top[2] == 1 && (i-=1)
  ra.residuals[1:i,:]
end


\(f::Function, b::Vector) = f(b)

#Borders are always white, issue on UnicodePlots? light terminals suffer
function showplot(vals::Vector)
  isdefined(Main, :UnicodePlots) || warn("UnicodePlots not found; no plotsies T.T ")
  println(lineplot(1:length(vals), vals, title = "Convergence", name = "resnorm"))
end

check(tol::Real, resnorm::Real, ::Vector, ::Integer, ::Type{Val{true}}) = resnorm < tol

function check(tol::Real, resnorm::Real, resnorms::Vector, iter::Int, ::Type{Val{false}})
  resnorms[iter] = resnorm
  if resnorms[iter] < tol
    resize!(resnorms,iter)
    return true
  end
  false
end

check(tol::Real, resnorm::Real, resnorms::Vector, iter::Integer) =
  check(tol,resnorm,resnorms,iter,Val{length(resnorms)==1})

#### Type-handling
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))
Amultype(A, x) = typeof(one(eltype(A))*one(eltype(x)))
if VERSION < v"0.4.0-dev+6068"
    real{T<:Real}(::Type{Complex{T}}) = T
    real{T<:Real}(::Type{T}) = T
end
eps{T<:Real}(::Type{Complex{T}}) = eps(T)

function randx(A, b)
    T = Adivtype(A, b)
    x = initrand!(Array(T, size(A, 2)))
end

function zerox(A, b)
    T = Adivtype(A, b)
    x = zeros(T, size(A, 2))
end

#### Numerics
function update!(x, α::Number, p::AbstractVector)
    for i = 1:length(x)
        x[i] += α*p[i]
    end
    x
end

function initrand!(v::Vector)
    _randn!(v)
    nv = norm(v)
    for i = 1:length(v)
        v[i] /= nv
    end
    v
end
_randn!(v::Array{Float64}) = randn!(v)
_randn!(v) = copy!(v, randn(length(v)))

#### Reporting
type ConvergenceHistory{T, R}
    isconverged::Bool
    threshold::T
    mvps::Int
    residuals::R
end

function empty!(ch::ConvergenceHistory)
    ch.isconverged = false
    ch.mvps = 0
    empty!(ch.residuals)
    ch
end

function push!(ch::ConvergenceHistory, resnorm::Number)
    push!(ch.residuals, resnorm)
    ch
end
push!(ch::ConvergenceHistory, residual::AbstractVector) = push!(ch, norm(residual))

#### Errors
export PosSemidefException

type PosSemidefException <: Exception
    msg :: AbstractString
    PosSemidefException(msg::AbstractString="Matrix was not positive semidefinite") = new(msg)
end

#### Support for linear-operators-defined-as-functions
#
# Suppose you have a function implementing multiplication of b by A.
# This function can have any syntax, but for the purposes of
# illustration let's suppose it's defined as
#   mulbyA!(output, b, Adata)
# Where Adata might be some parameters that your function needs.
# You can represent it as a linear operator using
#   A = MatrixFcn{T}(m, n, (output, b) -> mulbyA!(output, b, Adata))
# where T is the "element type" of Adata.
# Note that there are a couple of requirements:
#   - mulbyA! stores the result in the pre-allocated output
#   - mulbyA! should also return output as its sole return value

# If the algorithm also needs multiplication by A', use MatrixCFcn
# instead and initialize it with both functions.

export MatrixFcn, MatrixCFcn

abstract AbstractMatrixFcn{T}

type MatrixFcn{T} <: AbstractMatrixFcn{T}
    m::Int
    n::Int
    mul::Function
end

type MatrixCFcn{T} <: AbstractMatrixFcn{T}
    m::Int
    n::Int
    mul::Function
    mulc::Function
end

MatrixFcn(A::AbstractMatrix) = MatrixFcn{eltype(A)}(size(A, 1), size(A, 2), (output, x)-> A_mul_B!(output, A, x))
MatrixCFcn(A::AbstractMatrix) = MatrixFcn{eltype(A)}(size(A, 1), size(A, 2), (output, x)->A_mul_B!(output, A, x), (output, b) -> Ac_mul_B(output, A, b))

eltype{T}(op::AbstractMatrixFcn{T}) = T

ndims(op::AbstractMatrixFcn) = 2

size(op::AbstractMatrixFcn, dim::Integer) = (dim == 1) ? op.m :
                                            (dim == 2) ? op.n : 1
size(op::AbstractMatrixFcn) = (op.m, op.n)

length(op::AbstractMatrixFcn) = op.m*op.n

*(op::AbstractMatrixFcn, b) = A_mul_B(op, b)

A_mul_B{R,S}(op::AbstractMatrixFcn{R}, b::AbstractVector{S}) = A_mul_B!(Array(promote_type(R,S), op.m), op, b)

A_mul_B!(output, op::AbstractMatrixFcn, b) = op.mul(output, b)

Ac_mul_B{R,S}(op::MatrixCFcn{R}, b::AbstractVector{S}) = Ac_mul_B!(Array(promote_type(R,S), op.n), op, b)

Ac_mul_B!(output, op::AbstractMatrixFcn, b) = op.mulc(output, b)
