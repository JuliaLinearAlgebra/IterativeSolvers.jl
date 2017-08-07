using IterativeSolvers
using Base.Test

@testset "rsvd_fnkz" begin

srand(1234321)

@testset "Rank 2 matrix" begin
    A = reshape(1 : 64, 8, 8)
    S = rsvd_fnkz(A, 2, ϵ=1e-9)
    @test vecnorm(A - S[:U] * Diagonal(S[:S]) * S[:Vt]) ≤ 1e-9
end

@testset "Linearly dependent matrix" begin
    A = [collect(1 : 8) collect(1 : 8) collect(1 : 8) collect(1 : 8)]
    S = rsvd_fnkz(A, 4)
    @test vecnorm(A - S[:U] * Diagonal(S[:S]) * S[:Vt]) ≤ 1e-9
end

#SVD of a matrix of zeros should be empty
@testset "Zero matrix" begin
    S = rsvd_fnkz(zeros(10, 10), 1)
    @test S[:U] == zeros(10, 0)
    @test S[:S] == []
    @test S[:Vt] == zeros(0, 10)
end

end
