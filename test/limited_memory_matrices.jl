using IterativeSolvers; const IS = IterativeSolvers
using LinearAlgebra
using SparseArrays
using Test

@testset "LimitedMemoryMatrix" begin
    A = IS.LimitedMemoryMatrix{Float64, Matrix{Float64}}(fill(1.0, 4, 4), 4, 4)
    @test A[1, 1] == 1
    IS.hcat!(A, fill(2.0, 4))
    @test_throws ArgumentError A[1, 1]
    @test A[:, end] == fill(2.0, 4)

    A = IS.LimitedMemoryMatrix(fill(1.0, 4), 4)
    @test_throws ArgumentError A[3, 2]
    @test A[3, 1] == 1
    IS.hcat!(A, fill(2.0, 4))
    @test A[3, 1] == 1
    @test A[3, 2] == 2
    @test_throws ArgumentError A[3, 3]
    IS.hcat!(A, fill(3.0, 4))
    IS.hcat!(A, fill(4.0, 4))
    IS.hcat!(A, fill(5.0, 4))
    @test_throws ArgumentError A[3, 1]
    @test A[3, 2] == 2
    @test A[3, 3] == 3
    @test A[3, 4] == 4
    @test A[3, 5] == 5
end