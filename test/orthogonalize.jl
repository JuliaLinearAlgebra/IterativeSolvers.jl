using IterativeSolvers
using FactCheck

srand(1234321)

facts("Orthogonalization") do
    n = 10
    m = 3

    for T = [Complex64, Float64]
        context("$T") do
    
        # Create an orthonormal matrix V
        V, = qr(rand(T, n, m))

        # And a random vector to be orth. to V.
        w_original = rand(T, n)

        for method = [DGKS, ClassicalGramSchmidt, ModifiedGramSchmidt]

            context("$method") do
            # Projection size
            h = zeros(T, m)

            # Orthogonalize w in-place
            w = copy(w_original)
            nrm = orthogonalize_and_normalize!(V, w, h, method)

            # Normality
            @fact norm(w) --> roughly(one(real(T)))

            # Orthogonality
            @fact norm(V' * w) --> roughly(zero(real(T)), atol = 10 * eps(real(T)))

            # Denormalizing and adding the components in V should give back the original
            @fact nrm * w + V * h --> roughly(w_original)

            end
        end
        end
    end
end
