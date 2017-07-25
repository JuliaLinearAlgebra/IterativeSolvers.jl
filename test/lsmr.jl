using IterativeSolvers
using Base.Test

import Base: size, A_mul_B!, Ac_mul_B!, eltype, similar, scale!, copy!, fill!, length
import Base.LinAlg: norm

srand(1234321)

# Type used in SOL test
type Wrapper
    m::Int
    n::Int
end

size(op::Wrapper, dim::Integer) = (dim == 1) ? op.m : (dim == 2) ? op.n : 1
size(op::Wrapper) = op.m, op.n
eltype(op::Wrapper) = Int

function A_mul_B!(α, A::Wrapper, x, β, y)
    m, n = size(A)
    scale!(y, β)
    y[1] = y[1] + α * x[1]
    for i = 2:n
        y[i] = y[i] + i * α * x[i] + (i-1) * α * x[i-1]
    end
    for i = n+1:m
        y[i] = y[i]
    end
    return y
end

function Ac_mul_B!(α, A::Wrapper, x, β, y)
    m, n = size(A)
    mn = min(m, n)
    scale!(y, β)
    for i = 1 : mn - 1
        y[i] = y[i] + α * i * (x[i] + x[i + 1])
    end
    y[mn] = y[mn] + α * mn * x[mn]
    for i = m + 1 : n
        y[i] = y[i]
    end
    return y
end

# Type used in Dampenedtest
# solve (A'A + diag(v).^2 ) x = b
# using LSMR in the augmented space A' = [A ; diag(v)] b' = [b; zeros(size(A, 2)]
type DampenedVector{Ty, Tx}
    y::Ty
    x::Tx
end

eltype(a::DampenedVector) = promote_type(eltype(a.y), eltype(a.x))
norm(a::DampenedVector) = sqrt(norm(a.y)^2 + norm(a.x)^2)

function copy!{Ty, Tx}(a::DampenedVector{Ty, Tx}, b::DampenedVector{Ty, Tx})
    copy!(a.y, b.y)
    copy!(a.x, b.x)
    a
end

function fill!(a::DampenedVector, α::Number)
    fill!(a.y, α)
    fill!(a.x, α)
    a
end

scale!(α::Number, a::DampenedVector) = scale!(a, α)

function scale!(a::DampenedVector, α::Number)
    scale!(a.y, α)
    scale!(a.x, α)
    a
end

similar(a::DampenedVector, T) = DampenedVector(similar(a.y, T), similar(a.x, T))
length(a::DampenedVector) = length(a.y) + length(a.x)


type DampenedMatrix{TA, Tx}
    A::TA
    diagonal::Tx
end

eltype(A::DampenedMatrix) = promote_type(eltype(A.A), eltype(A.diagonal))

function size(A::DampenedMatrix, dim::Integer)
    m, n = size(A.A)
    l = length(A.diagonal)
    dim == 1 ? (m + l) :
    dim == 2 ? n : 1
end

function A_mul_B!{TA, Tx, Ty}(α::Number, mw::DampenedMatrix{TA, Tx}, a::Tx,
                β::Number, b::DampenedVector{Ty, Tx})
    if β != 1.
        if β == 0.
            fill!(b, 0.)
        else
            scale!(b, β)
        end
    end
    A_mul_B!(α, mw.A, a, 1.0, b.y)
    map!((z, x, y)-> z + α * x * y, b.x, b.x, a, mw.diagonal)
    return b
end

function Ac_mul_B!{TA, Tx, Ty}(α::Number, mw::DampenedMatrix{TA, Tx}, a::DampenedVector{Ty, Tx},
                β::Number, b::Tx)
    if β != 1.
        if β == 0.
            fill!(b, 0.)
        else
            scale!(b, β)
        end
    end
    Ac_mul_B!(α, mw.A, a.y, 1.0, b)
    map!((z, x, y)-> z + α * x * y, b, b, a.x, mw.diagonal)
    return b
end

@testset "LSMR" begin

    @testset "Small dense matrix" begin
        A = rand(10, 5)
        b = rand(10)
        x = lsmr(A, b)
        @test norm(x - A\b) ≤ √eps()
    end

    @testset "SOL test" for (m, n, damp) = ((10, 10, 0), (20, 10, 0), (20, 10, 0.1))
        # Test adapted from the BSD-licensed Matlab implementation at
        #    http://www.stanford.edu/group/SOL/software/lsqr.html
        #              Michael Saunders, Systems Optimization Laboratory,
        #              Dept of MS&E, Stanford University.
        #-----------------------------------------------------------------------

        # This is a simple example for testing  LSMR.
        # It uses the leading m*n submatrix from
        # A = [ 1
        #       1 2
        #         2 3
        #           3 4
        #             ...
        #               n ]
        # suitably padded by zeros.
        #
        # 11 Apr 1996: First version for distribution with lsqr.m.
        #              Michael Saunders, Dept of EESOR, Stanford University.

        A = Wrapper(m, n)
        xtrue = n : -1 : 1
        b = Vector{Float64}(m)
        A_mul_B!(1.0, A, xtrue, 0.0, b)
        x = lsmr(A, b, atol = 1e-7, btol = 1e-7, conlim = 1e10, maxiter = 10n)
        r = A_mul_B!(-1, A, x, 1, b)
        @test norm(r) ≤ 1e-4
    end

    @testset "Dampened test" for (m, n) = ((10, 10), (20, 10))
        # Test used to make sure A, b can be generic matrix / vector
        srand(1234)
        
        b = rand(m)
        A = rand(m, n)
        v = rand(n)
        Adampened = DampenedMatrix(A, v)
        bdampened = DampenedVector(b, zeros(n))
        x, ch = lsmr(Adampened, bdampened, log=true)
        @test norm((A'A + diagm(v) .^ 2)x - A'b) ≤ 1e-3
    end
end
