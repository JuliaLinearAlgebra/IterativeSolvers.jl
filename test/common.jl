module TestCommon

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

@testset "Vector{$T}, conjugated and unconjugated dot products" for T in (ComplexF32, ComplexF64)
    n = 100
    x = rand(T, n)
    y = rand(T, n)

    # Conjugated dot product
    @test IterativeSolvers._dot(x, Val(true)) ≈ x'x
    @test IterativeSolvers._dot(x, y, Val(true)) ≈ x'y
    @test IterativeSolvers._dot(x, Val(true)) ≈ IterativeSolvers._dot(x, x, Val(true))

    # Unonjugated dot product
    @test IterativeSolvers._dot(x, Val(false)) ≈ transpose(x) * x
    @test IterativeSolvers._dot(x, y, Val(false)) ≈ transpose(x) * y
    @test IterativeSolvers._dot(x, Val(false)) ≈ IterativeSolvers._dot(x, x, Val(false))
end

end

DocMeta.setdocmeta!(IterativeSolvers, :DocTestSetup, :(using IterativeSolvers); recursive=true)
doctest(IterativeSolvers, manual=true)

end # module
