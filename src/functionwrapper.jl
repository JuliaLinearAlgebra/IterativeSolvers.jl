# If you have a function something like
#   mulbyA!(output, Adata, b)
# which implements multiplication by A, you can "convert" it to an
# AbstractMatrix using
#   A = MatrixFcn{T}(m, n, (output, b) -> mulbyA!(output, Adata, b))
# where T is the type of Adata
# Note that:
#   - mulbyA! stores the result in the pre-allocated output
#   - mulbyA! should also return output

# If the algorithm also needs multiplication by A', use MatrixCFcn
# instead and supply both functions.

import Base: size, *, A_mul_B, A_mul_B!, Ac_mul_B, Ac_mul_B!

export MatrixFcn, MatrixCFcn

abstract AbstractMatrixFcn{T} <: AbstractMatrix{T}

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

size(A::AbstractMatrixFcn, dim::Integer) = (dim == 1) ? A.m : (dim == 2) ? A.n : 1
size(A::AbstractMatrixFcn) = (A.m, A.n)

*(A::AbstractMatrixFcn, b) = A_mul_B(A, b)

A_mul_B{R,S}(A::AbstractMatrixFcn{R}, b::AbstractVector{S}) = A_mul_B!(Array(promote_type(R,S), A.m), A, b)

A_mul_B!(output, A::AbstractMatrixFcn, b) = A.mul(output, b)

Ac_mul_B{R,S}(A::MatrixCFcn{R}, b::AbstractVector{S}) = Ac_mul_B!(Array(promote_type(R,S), A.n), A, b)

Ac_mul_B!(output, A::AbstractMatrixFcn, b) = A.mulc(output, b)
