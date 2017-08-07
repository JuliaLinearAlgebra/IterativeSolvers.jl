using IterativeSolvers
using Base.Test

srand(1234321)

@testset "Basic operations" begin

@testset "Adivtype" begin
    A = rand(1:10, 5, 5)
    b = rand(Float32, 5)
    @test IterativeSolvers.Adivtype(A, b) == Float32
end

@testset "Identity preconditioner" begin
    P = Identity()
    x = rand(10)
    y = zeros(x)

    # Should be a no-op
    @test P \ x == x
    @test A_ldiv_B!(P, copy(x)) == x
    @test A_ldiv_B!(y, P, copy(x)) == x
end

end
