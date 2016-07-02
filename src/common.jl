using UnicodePlots
import Base: \, eltype, empty!, eps, length, ndims, push!, real, size, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export A_mul_B

\(f::Function, b::Vector) = f(b)

showplot(vals::Residuals) = println(buildplot(vals))

function buildplot(ra::ResArray)
    vals = ra.residuals
    lineplot(collect(1:length(vals)), vals)
end
function buildplot(ra::RestResArray)
    vals = ra.residuals
    maxy = maximum(vals)
    miny = minimum(vals)
    iters = ra.restart
    restarts = div(length(vals),iters)
    plot = lineplot([0],[0],xlim=[1,restarts*iters],ylim=[miny,maxy])
    for restart in 1:restarts
        first = iters*(restart-1)
        last = restart*iters-1
        lineplot!(plot,collect(first:last),collect(vals[first:last]), color=:blue)
        lineplot!(plot,[first,first],[miny,maxy], color=:white)
    end
    plot
end

abstract Residuals{T,N}
typealias SolResidues Residuals{T,1}
typealias EigResidues Residuals{T,2}
type ResSingle{T,N} <: Residuals{T,N}
  residuals::{T,N}
end
type ResArray{T,N} <: Residuals{T,N}
  top::Int
  residuals::Array{T,N}
end
type RestResArray{T,N} <: Residuals{T,N}
  top::Int
  restart::Int
  residuals::Array{T,N}
end
ResSingle() = ResSingle(0)
ResArray(maxiter) = ResArray(1,zeros(maxiter))
RestResArray(maxiter,restart) = RestResArray(1,restart,zeros(maxiter*restart))
push!(rs::ResSingle, res::Real) = (rs.residual = res; nothing)
function push!(ra::Union{ResArray,RestResArray}, res::Real)
  ra.residuals[ra.top] = res
  ra.top+=1
  nothing
end
shrink!(rs::ResSingle) = ra.residual
shrink!(ra::Union{ResArray,RestResArray}) = resize!(ra.residuals,ra.top-1)

last(rs::ResSingle) = rs.residual
last(ra::Union{ResArray,RestResArray}) = ra.residuals[ra.top-1]

iters(r::Union{ResArray,RestResArray}) = ra.top-1
isconverged(r::Residuals,tol::Real) = 0 <= last(r) < tol

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
type ConvergenceHistory{T}
    isconverged::Bool
    threshold::T
    mvps::Int
    residuals::Residuals
end
ConvergenceHistory{T}(res::Residuals,tol::T,mvps::Int) =
    ConvergenceHistory(isconverged(res,tol),tol,mvps,res)

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
