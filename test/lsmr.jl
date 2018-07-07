using IterativeSolvers
using Test
using LinearAlgebra
using Random
using SparseArrays

import Base: size, eltype, similar, copyto!, fill!, length
import LinearAlgebra: norm, mul!, rmul!, lmul!

# Type used in Dampenedtest
# solve (A'A + diag(v).^2 ) x = b
# using LSMR in the augmented space A' = [A ; diag(v)] b' = [b; zeros(size(A, 2)]
mutable struct DampenedVector{Ty, Tx}
    y::Ty
    x::Tx
end

eltype(a::DampenedVector) = promote_type(eltype(a.y), eltype(a.x))
norm(a::DampenedVector) = sqrt(norm(a.y)^2 + norm(a.x)^2)

function Base.Broadcast.broadcast!(f::Tf, to::DampenedVector, from::DampenedVector, args...) where {Tf}
    to.x .= f.(from.x, args...)
    to.y .= f.(from.y, args...)
    to
end

function copyto!(a::DampenedVector{Ty, Tx}, b::DampenedVector{Ty, Tx}) where {Ty, Tx}
    copyto!(a.y, b.y)
    copyto!(a.x, b.x)
    a
end

function fill!(a::DampenedVector, α::Number)
    fill!(a.y, α)
    fill!(a.x, α)
    a
end

lmul!(α::Number, a::DampenedVector) = rmul!(a, α)

function rmul!(a::DampenedVector, α::Number)
    rmul!(a.y, α)
    rmul!(a.x, α)
    a
end

similar(a::DampenedVector, T) = DampenedVector(similar(a.y, T), similar(a.x, T))
length(a::DampenedVector) = length(a.y) + length(a.x)

mutable struct DampenedMatrix{TA, Tx}
    A::TA
    diagonal::Tx
end

eltype(A::DampenedMatrix) = promote_type(eltype(A.A), eltype(A.diagonal))

function size(A::DampenedMatrix)
    m, n = size(A.A)
    l = length(A.diagonal)
    (m + l, n)
end

function size(A::DampenedMatrix, dim::Integer)
    m, n = size(A.A)
    l = length(A.diagonal)
    dim == 1 ? (m + l) :
    dim == 2 ? n : 1
end

function mul!(b::DampenedVector{Ty, Tx}, mw::DampenedMatrix{TA, Tx}, a::Tx,
    α::Number, β::Number) where {TA, Tx, Ty}
    if β != 1.
        if β == 0.
            fill!(b, 0.)
        else
            rmul!(b, β)
        end
    end
    mul!(b.y, mw.A, a, α, 1.0)
    map!((z, x, y)-> z + α * x * y, b.x, b.x, a, mw.diagonal)
    return b
end

function mul!(b::Tx, mw::Adjoint{DampenedMatrix{TA, Tx}}, a::DampenedVector{Ty, Tx},
    α::Number, β::Number) where {TA, Tx, Ty}
    if β != 1.
        if β == 0.
            fill!(b, 0.)
        else
            rmul!(b, β)
        end
    end
    mul!(b, adjoint(mw.A), a.y, α, 1.0)
    map!((z, x, y)-> z + α * x * y, b, b, a.x, mw.diagonal)
    return b
end

"""
Produces the m × n submatrix from
A = [ 1
      1 2
        2 3
          3 4
            ...
              n ]
suitably padded by zeros.
"""
function sol_matrix(m, n)
    mn = min(m, n)
    spdiagm((1.0 : mn - 1, 1.0 : mn), (-1, 0), m, n)
end

@testset "LSMR" begin
    srand(1234321)

    @testset "Small dense matrix" for T = (Float32, Float64)
        A = rand(T, 10, 5)
        b = rand(T, 10)
        x, history = lsmr(A, b, log = true)
        @test isa(history, ConvergenceHistory)
        @test norm(x - A\b) ≤ √eps(T)
    end

    @testset "SOL test" for (m, n, damp) = ((10, 10, 0), (20, 10, 0), (20, 10, 0.1))
        # Test adapted from the BSD-licensed Matlab implementation at
        #    http://www.stanford.edu/group/SOL/software/lsqr.html
        #              Michael Saunders, Systems Optimization Laboratory,
        #              Dept of MS&E, Stanford University.
        #-----------------------------------------------------------------------
        # 11 Apr 1996: First version for distribution with lsqr.m.
        #              Michael Saunders, Dept of EESOR, Stanford University.

        A = sol_matrix(m, n)
        x = float(n : -1 : 1)
        b = A * x
        x_lsmr = lsmr(A, b, atol = 1e-7, btol = 1e-7, conlim = 1e10, maxiter = 10n)
        @test norm(b - A * x) ≤ 1e-4
    end

    # @testset "Dampened test" for (m, n) = ((10, 10), (20, 10))
    #     # Test used to make sure A, b can be generic matrix / vector
    #     b = rand(m)
    #     A = rand(m, n)
    #     v = rand(n)
    #     Adampened = DampenedMatrix(A, v)
    #     bdampened = DampenedVector(b, zeros(n))
    #     x, ch = lsmr(Adampened, bdampened, log=true)
    #     @test norm((A'A + Matrix(Diagonal(v)) .^ 2)x - A'b) ≤ 1e-3
    # end
end
