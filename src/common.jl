import Base: eltype, empty!, length, ndims, push!, size, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export A_mul_B

#### Type-handling
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))
Amultype(A, x) = typeof(one(eltype(A))*one(eltype(x)))

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

########################
# Termination criteria #
########################

@doc doc"""
Abstract termination criterion
""" ->
abstract stopcriterion

@doc doc"""
Termination criterion state

Fields:

    iter: Iteration number
    resnorm²: square of the residual norm
""" ->
immutable stopcriteriastate{T}
    iter :: Int
    resnorm² :: T
end



@doc doc"""
Stopping criterion based on maximum number of iterations

Field:

    maxiter: maximum number of iterations allowed
""" ->
immutable maxiterations <: stopcriterion
    maxiter :: Int
end

@doc doc"""
Check if maximum number of iterations has been reached or exceeded.
"""
done(criterion::maxiterations, current::stopcriteriastate) = 
    current.iter ≥ criterion.maxiter



@doc doc"""
Stopping criterion based on absolute value of the residual norm

Field:

    threshold: value of threshold
""" ->
immutable absresnorm{T<:Real} <: stopcriterion
    threshold :: T
end

@doc doc"""
Check if the residual norm has fallen below a certain absolute threshold.
""" ->
done(criterion::absresnorm, current::stopcriteriastate) =
    current.resnorm²<criterion.threshold^2



@doc doc"""
Stopping criterion based on relative value of the residual norm

Field:

    relthreshold: value of relative threshold

Implementation note:

    A relresnorm stopcriterion gets replaced by an absresnorm
    once a stopcriterionstate is computed. In other words, the
    continuation of a relresnorm is an absresnorm once the initial resnorm is
    computed.
""" ->
immutable relresnorm{T<:Real} <: stopcriterion
    relthreshold :: T
end

@doc doc"""
Check if the residual norm has fallen below a certain relative threshold.

Implementation note:

    Always return false. A relresnorm stopcriterion gets replaced by an absnorm
    stopcriterion once the initial resnorm is known. Hence, if this criterion
    is still around, the starting norm is unknown and there's no way to know if
    it should stop.
""" ->
done(criterion::relresnorm, current::stopcriteriastate) =
    false



@doc doc"""
A type encapsulating various termination criteria.
""" ->
type Terminator
    criteria :: Vector{stopcriterion}
end

@doc doc"""
Add a stopcriterion.
""" ->
push!(T::Terminator, criterion::stopcriterion) =
    push!(T.criteria, criterion)

@doc doc"""
    Check whether to terminate algorithm based on at least one of the
    stopcriteria returning `true`.
""" ->
function done(termination::Terminator, currentstate::stopcriteriastate)
    for (idx, criterion) in enumerate(termination.criteria)
        if isa(criterion, relresnorm) && isfinite(currentstate.resnorm²)
            #There is a computed resnorm, so replace the
            #relative criterion it by its continuation as an absresnorm
            #Be careful to preserve ordering of criteria
            #since we are modifying the list while iterating over it
            deleteat!(termination.criteria, idx)
            insert!(termination.criteria, idx,
                absresnorm(criterion.relthreshold*√currentstate.resnorm²))
        else #check if current criterion is good to terminate
            done(criterion, currentstate) && return true
        end
    end
    return length(termination.criteria)==0 #If no criteria, always terminate
end



# Algorithms

@doc doc"""
Abstract iterative solver
""" ->
abstract IterativeSolver

@doc doc"""
Abstract state of iterative solver
""" ->
abstract IterationState

immutable KrylovSpace #XXX to replace KrylovSubspace
    A
    v0 :: Vector
end

