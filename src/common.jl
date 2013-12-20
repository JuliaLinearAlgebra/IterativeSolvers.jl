import Base: eltype, empty!, length, ndims, push!, size, *, A_mul_B, A_mul_B!, Ac_mul_B, Ac_mul_B!

#### Type-handling
Abtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))

#### Reporting
type ConvergenceHistory{T, R}
    isconverged::Bool
    threshold::T
    mvps::Int
    residuals::R
end

ConvergenceHistory{T}(isconverged::Bool, threshold::T, mvps::Int, resnorms::AbstractVector{T}) = ConvergenceHistory{T,Vector{T}}(isconverged, threshold, mvps, convert(Vector, resnorms))
ConvergenceHistory{T}(isconverged::Bool, threshold, mvps::Integer, resnorms::AbstractVector{T}) = ConvergenceHistory{T,Vector{T}}(isconverged, convert(T, threshold), int(mvps), convert(Vector, resnorms))
ConvergenceHistory{T}(::Type{T}) = ConvergenceHistory{T,Vector{T}}(false, zero(T), 0, T[])

function empty!(ch::ConvergenceHistory)
    ch.isconverged = false
    ch.threshold = 0
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
    msg :: UTF8String
    PosSemidefException(msg::String="Matrix was not positive semidefinite") = new(msg)
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
