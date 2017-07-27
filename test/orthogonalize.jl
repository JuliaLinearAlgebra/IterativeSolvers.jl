using IterativeSolvers
using Base.Test

@testset "Orthogonalization" begin

srand(1234321)
n = 10
m = 3

@testset "Eltype $T" for T = (Complex64, Float64)

    # Create an orthonormal matrix V
    V, = qr(rand(T, n, m))

    # And a random vector to be orth. to V.
    w_original = rand(T, n)

    @testset "Using $method" for method = (DGKS, ClassicalGramSchmidt, ModifiedGramSchmidt)

        # Projection size
        h = zeros(T, m)

        # Orthogonalize w in-place
        w = copy(w_original)
        nrm = orthogonalize_and_normalize!(V, w, h, method)

        # Normality
        @test norm(w) ≈ one(real(T))

        # Orthogonality
        @test norm(V'w) ≈ zero(real(T)) atol = 10eps(real(T))

        # Denormalizing and adding the components in V should give back the original
        @test nrm * w + V * h ≈ w_original
    end
end
end
