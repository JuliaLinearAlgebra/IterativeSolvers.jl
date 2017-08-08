import Base.LinAlg: A_mul_B!
import Base: start, next, done

struct OffDiagonal{matT <: SparseMatrixCSC}
    data::matT
end

struct StrictlyUpperTriangular{matT <: SparseMatrixCSC}
    data::matT
end

"""
FastLowerTriangular is a data structure that
pre-computes the indices of the diagonal elements
so that backward substitution can be performed without searching.
Throws when the lower-triangular matrix is nonsingular.
"""
struct FastLowerTriangular{matT <: SparseMatrixCSC}
    data::matT
    diag_idx::Vector{Int}

    function FastLowerTriangular{matT}(A::matT) where {matT}
        # Check square?
        diag_idx = Vector{Int}(A.n)

        for col = 1 : A.n
            r1 = A.colptr[col]
            r2 = A.colptr[col + 1] - 1
            r1 = searchsortedfirst(A.rowval, col, r1, r2, Base.Order.Forward)
            if (r1 > r2) || (A.rowval[r1] != col)
                throw(Base.LinAlg.SingularException(col))
            end
            diag_idx[col] = r1
        end

        new(A, diag_idx)
    end
end

FastLowerTriangular(A::matT) where {matT<:SparseMatrixCSC} = FastLowerTriangular{matT}(A)

"""
Backward substitution for the FastLowerTriangular type
"""
function back_sub!(F::FastLowerTriangular, x::StridedVector{T}) where {T}
    A = F.data

    @inbounds for col = 1 : A.n

        # Solve for diagonal element
        idx = F.diag_idx[col]
        x[col] = x[col] / A.nzval[idx]

        # Substitute next values involving x[col]
        for i = idx + 1 : (A.colptr[col + 1] - 1)
            x[A.rowval[i]] -= A.nzval[i] * x[col]
        end
    end
end

"""
Backward substitution with α = 1 and β = 0. The α, β and y parameters
are useful for SOR.
"""
function back_sub!(α, F::FastLowerTriangular, x::StridedVector, β, y::StridedVector)
    A = F.data

    @inbounds for col = 1 : A.n

        # Solve for diagonal element
        idx = F.diag_idx[col]
        x[col] = α * x[col] / A.nzval[idx] + β * y[col]

        # Substitute next values involving x[col]
        for i = idx + 1 : (A.colptr[col + 1] - 1)
            x[A.rowval[i]] -= A.nzval[i] * x[col]
        end
    end
end

function A_mul_B!(α, A::OffDiagonal, x::StridedVector, β, y::StridedVector)
    # Specialize for β = 0 and β = 1
    M = A.data

    if β != 1
        if iszero(β)
            fill!(y, zero(eltype(y)))
        else
            scale!(β, y)
        end
    end

    @inbounds for col = 1 : M.n
        αx = α * x[col]
        for j = M.colptr[col] : (M.colptr[col + 1] - 1)
            # Skip the diagonal.
            if col != M.rowval[j]
                y[M.rowval[j]] += M.nzval[j] * αx
            end
        end
    end

    y
end

"""
Computes z := α * A * x + β * y. Because A is StrictlyUpperTriangular
one can set z = x and update x in-place as x := α * A * x + β * y.
"""
function gauss_seidel_multiply!(α::T, A::StrictlyUpperTriangular, x::StridedVector{T}, β::T, y::StridedVector{T}, z::StridedVector{T}) where {T}
    M = A.data

    for col = 1 : M.n
        αx = α * x[col]
        @inbounds for j = M.colptr[col] : (M.colptr[col + 1] - 1)
            # Skip off-diagonal
            if M.rowval[j] ≥ col
                break
            end

            z[M.rowval[j]] += M.nzval[j] * αx
        end
        z[col] = β * y[col]
    end
    z
end

mutable struct JacobiIterable{vecT <: AbstractVector}
    R::OffDiagonal
    D::Diagonal

    x::vecT
    next::vecT
    b::vecT

    maxiter::Int
end

start(::JacobiIterable) = 1
done(j::JacobiIterable, iteration::Int) = iteration > j.maxiter
function next(j::JacobiIterable, iteration::Int)
    # tmp = D \ (b - R * x)
    T = eltype(j.x)
    copy!(j.next, j.b)
    A_mul_B!(-one(T), j.R, j.x, one(T), j.next)
    A_ldiv_B!(j.D, j.next)

    # Swap
    j.x, j.next = j.next, j.x

    nothing, iteration + 1
end

jacobi_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10) =
    JacobiIterable(OffDiagonal(A), Diagonal(A), x, similar(x), b, maxiter)

"""
Jacobi iteration is a matrix-splitting iteration where
`A = D + R` where `D` is diagonal and `R` the off-diagonal entries.
It iterates as x ← inv(D) * (b - R * x).
Storage costs are O(2n): a temporary vector is used and the diagonal
is constructed prior to the iterative process.
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

gauss_seidel_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10) = 
    GaussSeidelIterable(StrictlyUpperTriangular(A), FastLowerTriangular(A), x, b, maxiter)

start(::GaussSeidelIterable) = 1
done(g::GaussSeidelIterable, iteration::Int) = iteration > g.maxiter
function next(g::GaussSeidelIterable, iteration::Int)
    # x ← L \ (-U * x + b)
    T = eltype(g.x)
    gauss_seidel_multiply!(-one(T), g.U, g.x, one(T), g.b, g.x)
    back_sub!(g.L, g.x)

    nothing, iteration + 1
end

"""
Gauss-Seidel iteration is a matrix-splitting iteration where
`A = L + U` with `L` lower-triangular and `U` strictly upper-triangular.
It iterates as x ← inv(L) * (b - U * x).
One additional vector of 
"""
function gauss_seidel!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector; maxiter::Int = 10)
    iterable = gauss_seidel_iterable(x, A, b, maxiter = maxiter)
    for item = iterable end
    iterable.x
end

##
## SOR
##

mutable struct SORIterable{vecT <: AbstractVector, numT <: Real}
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
function next(s::SORIterable, iteration::Int)
    T = eltype(s.x)

    # next = b - U * x
    gauss_seidel_multiply!(-one(T), s.U, s.x, one(T), s.b, s.next)

    # next = ω * inv(L) * next + (1 - ω) * x
    back_sub!(s.ω, s.L, s.next, one(T) - s.ω, s.x)

    # Switch current and next iterate
    s.x, s.next = s.next, s.x

    nothing, iteration + 1
end

sor_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10) = 
    SORIterable(StrictlyUpperTriangular(A), FastLowerTriangular(A), ω, x, similar(x), b, maxiter)

"""
SOR iteration uses weights to control overshooting in Gauss-Seidel:
```
for i = 1 : n
  xᵢ_new ← (1 - ω)xᵢ_old + ω gauss_seidel_iterate_xᵢ
end
```
Two additional vectors are stored: a 
"""
function sor!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10)
    iterable = sor_iterable(x, A, b, ω, maxiter = maxiter)
    for item = iterable end
    iterable.x
end


mutable struct SSORIterable{matT <: SparseMatrixCSC, vecT, numT <: Real}
    A::matT
    ω::numT
    x::vecT
    b::vecT
    n::Int
    maxiter::Int
end

ssor_iterable(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10) = 
    SSORIterable(A, ω, x, b, length(x), maxiter)

start(s::SSORIterable) = 1
done(s::SSORIterable, iteration::Int) = iteration > s.maxiter

function next(s::SSORIterable, iteration::Int)
    A, x, b = s.A, s.x, s.b

    @inbounds for col = 1 : s.n
        diag_el = σ = zero(eltype(A))

        for idx = A.colptr[col] : A.colptr[col + 1] - 1
            if A.rowval[idx] == col
                diag_el = A.nzval[idx]
            else
                σ += A.nzval[idx] * x[A.rowval[idx]]
            end
        end

        x[col] += s.ω * ((b[col] - σ) / diag_el - x[col])
    end

    @inbounds for col = s.n : -1 : 1
        diag_el = σ = zero(eltype(A))

        for idx = A.colptr[col] : A.colptr[col + 1] - 1
            if A.rowval[idx] == col
                diag_el = A.nzval[idx]
            else
                σ += A.nzval[idx] * x[A.rowval[idx]]
            end
        end

        x[col] += s.ω * ((b[col] - σ) / diag_el - x[col])
    end

    nothing, iteration + 1
end

"""
Symmetric SOR. Assumes `A` is symmetric with all its
elements stored.

All operations are in-place in x.
"""
function ssor!(x::AbstractVector, A::SparseMatrixCSC, b::AbstractVector, ω::Real; maxiter::Int = 10)
    iterable = ssor_iterable(x, A, b, ω, maxiter = maxiter)
    for item = iterable end
    iterable.x
end