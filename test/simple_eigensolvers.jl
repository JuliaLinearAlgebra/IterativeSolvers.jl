using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

#######################
# Eigensystem solvers #
#######################

facts("Simple eigensolvers") do

n = 10

for T in (Float32, Float64, Complex64, Complex128)

    context("Matrix{$T}") do

    A = rand(T, n, n)
    A = A' * A
    λs = eigvals(A)

    tol = (eltype(T) <: Complex ? 2 : 1) * n^2 * cond(A) * eps(real(one(T)))

    ## Simple methods

    context("Power iteration") do
        λ, x = powm(A; tol = tol, maxiter = 10n)
        @fact λs[end] --> roughly(λ)
        @fact norm(A * x - λ * x) --> less_than(tol)
    end

    context("Inverse iteration") do
        # Set a target near the middle eigenvalue
        idx = div(n, 2)
        σ = T(0.75 * λs[idx] + 0.25 * λs[idx + 1])

        # Construct F = inv(A - σI) "matrix free"
        # Make sure we use complex arithmetic everywhere,
        # because of the follow bug in base: https://github.com/JuliaLang/julia/issues/22683
        F = lufact(complex(A) - UniformScaling(σ))
        Fmap = LinearMap((y, x) -> A_ldiv_B!(y, F, x), size(A, 1), complex(T), ismutating = true)

        λ, x = invpowm(Fmap; shift = σ, tol = tol, maxiter = 10n)

        @fact norm(A * x - λ * x) --> less_than(tol)
        @fact λ --> roughly(λs[idx])

    end

    end
end
end
