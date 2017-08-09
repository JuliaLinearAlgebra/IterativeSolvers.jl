import Base.LinAlg: A_mul_B!, A_ldiv_B!
import Base: start, next, done, getindex

struct DiagonalIndices{Tv, Ti <: Integer}
    matrix::SparseMatrixCSC{Tv,Ti}
    diag::Vector{Ti}

    function DiagonalIndices{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
        # Check square?
        diag = Vector{Ti}(A.n)

        for col = 1 : A.n
            r1 = A.colptr[col]
            r2 = A.colptr[col + 1] - 1
            r1 = searchsortedfirst(A.rowval, col, r1, r2, Base.Order.Forward)
            if r1 > r2 || A.rowval[r1] != col || iszero(A.nzval[r1])
                throw(Base.LinAlg.SingularException(col))
            end
            diag[col] = r1
        end

        new(A, diag)
    end
end

DiagonalIndices(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = DiagonalIndices{Tv,Ti}(A)

function A_ldiv_B!(y::AbstractVector{Tv}, D::DiagonalIndices{Tv,Ti}, x::AbstractVector{Tv}) where {Tv,Ti}
    @inbounds for row = 1 : D.matrix.n
        y[row] = x[row] / D.matrix.nzval[D.diag[row]]
    end
    y
end

@inline getindex(d::DiagonalIndices, i::Int) = d.diag[i]

struct FastLowerTriangular{Tv,Ti}
    matrix::SparseMatrixCSC{Tv,Ti}
    diag::DiagonalIndices{Tv,Ti}
end

struct FastUpperTriangular{Tv,Ti}
    matrix::SparseMatrixCSC{Tv,Ti}
    diag::DiagonalIndices{Tv,Ti}
end

struct StrictlyUpperTriangular{Tv,Ti}
    matrix::SparseMatrixCSC{Tv,Ti}
    diag::DiagonalIndices{Tv,Ti}
end

struct StrictlyLowerTriangular{Tv,Ti}
    matrix::SparseMatrixCSC{Tv,Ti}
    diag::DiagonalIndices{Tv,Ti}
end

struct OffDiagonal{Tv,Ti}
    matrix::SparseMatrixCSC{Tv,Ti}
    diag::DiagonalIndices{Tv,Ti}
end

"""
Forward substitution for the FastLowerTriangular type
"""
function forward_sub!(F::FastLowerTriangular, x::AbstractVector)
    A = F.matrix

    @inbounds for col = 1 : A.n

        # Solve for diagonal element
        idx = F.diag[col]
        x[col] = x[col] / A.nzval[idx]

        # Substitute next values involving x[col]
        for i = idx + 1 : (A.colptr[col + 1] - 1)
            x[A.rowval[i]] -= A.nzval[i] * x[col]
        end
    end

    x
end

"""
Forward substitution
"""
function forward_sub!(α, F::FastLowerTriangular, x::AbstractVector, β, y::AbstractVector)
    A = F.matrix

    @inbounds for col = 1 : A.n

        # Solve for diagonal element
        idx = F.diag[col]
        x[col] = α * x[col] / A.nzval[idx] + β * y[col]

        # Substitute next values involving x[col]
        for i = idx + 1 : (A.colptr[col + 1] - 1)
            x[A.rowval[i]] -= A.nzval[i] * x[col]
        end
    end

    x
end

"""
Backward substitution for the FastUpperTriangular type
"""
function backward_sub!(F::FastUpperTriangular, x::AbstractVector)
    A = F.matrix

    @inbounds for col = A.n : -1 : 1

        # Solve for diagonal element
        idx = F.diag[col]
        x[col] = x[col] / A.nzval[idx]

        # Substitute next values involving x[col]
        for i = A.colptr[col] : idx - 1
            x[A.rowval[i]] -= A.nzval[i] * x[col]
        end
    end

    x
end

function backward_sub!(α, F::FastUpperTriangular, x::AbstractVector, β, y::AbstractVector)
    A = F.matrix

    @inbounds for col = A.n : -1 : 1

        # Solve for diagonal element
        idx = F.diag[col]
        x[col] = α * x[col] / A.nzval[idx] + β * y[col]

        # Substitute next values involving x[col]
        for i = A.colptr[col] : idx - 1
            x[A.rowval[i]] -= A.nzval[i] * x[col]
        end
    end

    x
end

"""
Do A_mul_B! with the off-diagonal elements of a matrix.
"""
function A_mul_B!(α::T, O::OffDiagonal, x::AbstractVector, β::T, y::AbstractVector) where {T}
    # Specialize for β = 0 and β = 1
    A = O.matrix

    if β != one(T)
        if iszero(β)
            fill!(y, zero(T))
        else
            scale!(β, y)
        end
    end

    @inbounds for col = 1 : A.n
        αx = α * x[col]
        diag_index = O.diag[col]
        for j = A.colptr[col] : diag_index - 1
            y[A.rowval[j]] += A.nzval[j] * αx
        end
        for j = diag_index + 1 : A.colptr[col + 1] - 1
            y[A.rowval[j]] += A.nzval[j] * αx
        end
    end

    y
end

"""
Computes z := α * U * x + β * y. Because U is StrictlyUpperTriangular
one can set z = x and update x in-place as x := α * U * x + β * y.
"""
function gauss_seidel_multiply!(α, U::StrictlyUpperTriangular, x::AbstractVector, β, y::AbstractVector, z::AbstractVector)
    A = U.matrix

    for col = 1 : A.n
        αx = α * x[col]
        diag_index = U.diag[col]
        @inbounds for j = A.colptr[col] : diag_index - 1
            z[A.rowval[j]] += A.nzval[j] * αx
        end
        z[col] = β * y[col]
    end
    z
end

"""
Computes z := α * L * x + β * y. Because A is StrictlyLowerTriangular
one can set z = x and update x in-place as x := α * L * x + β * y.
"""
function gauss_seidel_multiply!(α, L::StrictlyLowerTriangular, x::AbstractVector, β, y::AbstractVector, z::AbstractVector)
    A = L.matrix

    for col = A.n : -1 : 1
        αx = α * x[col]
        z[col] = β * y[col]
        @inbounds for j = L.diag[col] + 1 : (A.colptr[col + 1] - 1)
            z[A.rowval[j]] += A.nzval[j] * αx
        end
    end
    z
end

##
## Jacobi
##

mutable struct JacobiIterable{T <: Number, vecT <: AbstractVector}
    O::OffDiagonal

    x::vecT
    next::vecT
    b::vecT

    maxiter::Int
end

start(::JacobiIterable) = 1
done(j::JacobiIterable, iteration::Int) = iteration > j.maxiter
function next(j::JacobiIterable{T}, iteration::Int) where {T}
    # tmp = D \ (b - (A - D) * x)
    copy!(j.next, j.b)
    A_mul_B!(-one(T), j.O, j.x, one(T), j.next)
    A_ldiv_B!(j.x, j.O.diag, j.next)

    nothing, iteration + 1
end

jacobi_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10) =
    JacobiIterable{eltype(x), typeof(x)}(OffDiagonal(A, DiagonalIndices(A)), x, similar(x), b, maxiter)

"""
Jacobi iteration is a matrix-splitting iteration with
`A = D + R` where `D` is the diagonal and `R` the off-diagonal.
It iterates as x ← D \\ (b - R * x).

The SparseMatrixCSC version of Jacobi iteration pre-computes the
indices of the diagonal of `A` and stores a temporary vector.
Storage costs are hence O(2n). For multiple solves with `A` it is
recommended to reuse the underlying iterator, see the 
`JacobiIterable` type.
"""
function jacobi!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10)
    iterable = jacobi_iterable(x, A, b, maxiter = maxiter)
    for item = iterable end
    iterable.x
end

##
## Gauss-Seidel
##

mutable struct GaussSeidelIterable{vecT <: AbstractVector}
    U::StrictlyUpperTriangular
    L::FastLowerTriangular

    x::vecT
    b::vecT

    maxiter::Int
end

function gauss_seidel_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10)
    D = DiagonalIndices(A)
    GaussSeidelIterable(StrictlyUpperTriangular(A, D), FastLowerTriangular(A, D), x, b, maxiter)
end

start(::GaussSeidelIterable) = 1
done(g::GaussSeidelIterable, iteration::Int) = iteration > g.maxiter
function next(g::GaussSeidelIterable, iteration::Int)
    # x ← L \ (-U * x + b)
    T = eltype(g.x)
    gauss_seidel_multiply!(-one(T), g.U, g.x, one(T), g.b, g.x)
    forward_sub!(g.L, g.x)

    nothing, iteration + 1
end

"""
Gauss-Seidel iteration is a matrix-splitting iteration with
`A = L + U` where `L` is lower-triangular and `U` strictly upper-triangular.
It iterates as x ← L \\ (b - U * x).

Gauss-Seidel can be computed in-place even for the `SparseMatrixCSC` type,
however, this implementation pre-computes the indices of the diagonal of `A`.
Hence storage costs are O(n). For multiple solves with `A` it is
recommended to reuse the underlying iterator, see the 
`GaussSeidelIterable` type.
"""
function gauss_seidel!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10)
    iterable = gauss_seidel_iterable(x, A, b, maxiter = maxiter)
    for item = iterable end
    iterable.x
end

##
## SOR
##

mutable struct SORIterable{T, vecT <: AbstractVector, numT <: Real}
    U::StrictlyUpperTriangular
    L::FastLowerTriangular
    ω::numT

    x::vecT
    next::vecT
    b::vecT

    maxiter::Int
end

start(::SORIterable) = 1
done(s::SORIterable, iteration::Int) = iteration > s.maxiter
function next(s::SORIterable{T}, iteration::Int) where {T}
    # next = b - U * x
    gauss_seidel_multiply!(-one(T), s.U, s.x, one(T), s.b, s.next)

    # next = ω * inv(L) * next + (1 - ω) * x
    forward_sub!(s.ω, s.L, s.next, one(T) - s.ω, s.x)

    # Switch current and next iterate
    s.x, s.next = s.next, s.x

    nothing, iteration + 1
end

function sor_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10)
    D = DiagonalIndices(A)
    T = eltype(x)
    SORIterable{T,typeof(x),eltype(ω)}(StrictlyUpperTriangular(A, D), FastLowerTriangular(A, D), ω, x, similar(x), b, maxiter)
end

"""
SOR takes a weighted average of the previous iterate and the Gauss-Seidel iterate
x_i ← ω * x_gs + (1 - ω) * x_i for each successive component i = 1, ..., n.

Although SOR can be computed fully in-place, this is difficult to accomplish for
the `SparseMatrixCSC` type. A temporary vector is pre-allocated together with the 
indices of the diagonal of `A`. Hence storage costs are O(2n). For multiple solves 
with `A` it is recommended to reuse the underlying iterator, see the 
`SORIterable` type.
"""
function sor!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10)
    iterable = sor_iterable(x, A, b, ω, maxiter = maxiter)
    for item = iterable end
    iterable.x
end

##
## SSOR
##

mutable struct SSORIterable{T, vecT, numT <: Real}
    sL::StrictlyLowerTriangular
    sU::StrictlyUpperTriangular
    L::FastLowerTriangular
    U::FastUpperTriangular
    ω::numT
    x::vecT
    tmp::vecT
    b::vecT
    maxiter::Int
end

function ssor_iterable(x::AbstractVector{T}, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10) where {T}
    D = DiagonalIndices(A)
    sL = StrictlyLowerTriangular(A, D)
    sU = StrictlyUpperTriangular(A, D)
    L = FastLowerTriangular(A, D)
    U = FastUpperTriangular(A, D)
    SSORIterable{T,typeof(x),typeof(ω)}(sL, sU, L, U, ω, x, similar(x), b, maxiter)
end

start(s::SSORIterable) = 1
done(s::SSORIterable, iteration::Int) = iteration > s.maxiter

function next(s::SSORIterable{T}, iteration::Int) where {T}
    # tmp = b - U * x
    gauss_seidel_multiply!(-one(T), s.sU, s.x, one(T), s.b, s.tmp)

    # tmp = ω * inv(L) * tmp + (1 - ω) * x
    forward_sub!(s.ω, s.L, s.tmp, one(T) - s.ω, s.x)

    # x = b - L * tmp
    gauss_seidel_multiply!(-one(T), s.sL, s.tmp, one(T), s.b, s.x)

    # x = ω * inv(U) * x + (1 - ω) * tmp
    backward_sub!(s.ω, s.U, s.x, one(T) - s.ω, s.tmp)

    nothing, iteration + 1
end

"""
Symmetric SOR (SSOR) combines two SOR sweeps (one forward and one backward) such that
the resulting iteration matrix is symmetric.

Although SSOR can be computed fully in-place, this is difficult to accomplish for
the `SparseMatrixCSC` type. A temporary vector is pre-allocated together with the 
indices of the diagonal of `A`. Hence storage costs are O(2n). For multiple solves 
with `A` it is recommended to reuse the underlying iterator, see the 
`SSORIterable` type.
"""
function ssor!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10)
    iterable = ssor_iterable(x, A, b, ω, maxiter = maxiter)
    for item = iterable end
    iterable.x
end
