using IterativeSolvers; const IS = IterativeSolvers
using LinearAlgebra
using SparseArrays
using Test

@testset "Limited Memory Matrices" begin

@testset "LimitedMemoryMatrix" begin
    A = IS.LimitedMemoryMatrix{Float64, Matrix{Float64}}(fill(1.0, 4, 4), 4, 4)
    @test A[1, 1] == 1
    IS.hcat!(A, fill(2.0, 4))
    @test_throws ArgumentError A[1, 1]
    @test A[:, end] == fill(2.0, 4)

    A = IS.LimitedMemoryMatrix(fill(1.0, 4), 4)
    @test A[3, 1] == 1
    IS.hcat!(A, fill(2.0, 4))
    @test A[3, 1] == 1
    @test A[3, 2] == 2
    IS.hcat!(A, fill(3.0, 4))
    IS.hcat!(A, fill(4.0, 4))
    IS.hcat!(A, fill(5.0, 4))
    @test_throws ArgumentError A[3, 1]
    @test A[3, 2] == 2
    @test A[3, 3] == 3
    @test A[3, 4] == 4
    @test A[3, 5] == 5
end

@testset "LimitedMemoryUpperTriangular" begin
    A = IS.LimitedMemoryUpperTriangular{Float64, Matrix{Float64}}(3)
    IS._grow_hcat!(A, fill(1.0, 1))
    IS._grow_hcat!(A, fill(2.0, 2))
    @test A ≈ [
        1.0 2.0
        0.0 2.0
    ]
    IS._grow_hcat!(A, fill(3.0, 3))
    v = fill(4.0, 4)
    v[1] = 0
    IS._grow_hcat!(A, v)
    @test A ≈ [
        0.0 2.0 3.0 0.0
        0.0 2.0 3.0 4.0
        0.0 0.0 3.0 4.0
        0.0 0.0 0.0 4.0
    ]
end

@testset "LimitedMemoryUpperHessenberg" begin
    A = IS.LimitedMemoryUpperHessenberg{Float64, Matrix{Float64}}(3)
    IS._grow_hcat!(A, fill(1.0, 2))
    IS._grow_hcat!(A, fill(2.0, 3))
    @test A ≈ [
        1.0 2.0    
        1.0 2.0
        0.0 2.0
    ]
    v = fill(3.0, 4)
    v[1] = 0
    IS._grow_hcat!(A, v)
    v = fill(4.0, 5)
    v[1:2] .= 0
    IS._grow_hcat!(A, v)
    @test A ≈ [
        0.0 2.0 0.0 0.0
        0.0 2.0 3.0 0.0
        0.0 2.0 3.0 4.0
        0.0 0.0 3.0 4.0
        0.0 0.0 0.0 4.0
    ]
end

end