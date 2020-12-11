using Documenter
using LinearAlgebra
using IterativeSolvers
using Test
using Random

Random.seed!(1234321)

@testset "Basic operations" begin

@testset "Adivtype" begin
    A = rand(1:10, 5, 5)
    b = rand(Float32, 5)
    @test IterativeSolvers.Adivtype(A, b) == Float32
end

@testset "Identity preconditioner" begin
    P = Identity()
    x = rand(10)
    y = zero(x)

    # Should be a no-op
    @test P \ x == x
    @test ldiv!(P, copy(x)) == x
    @test ldiv!(y, P, copy(x)) == x
end

end

DocMeta.setdocmeta!(IterativeSolvers, :DocTestSetup, :(using IterativeSolvers); recursive=true)
doctest(IterativeSolvers, manual=false)
