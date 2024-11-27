module TestLOBPCG

using IterativeSolvers
using LinearMaps
using LinearAlgebra
using Test
using Random
using SparseArrays


include("laplace_matrix.jl")

struct JacobiPrec{TD}
    diagonal::TD
end

LinearAlgebra.ldiv!(y, P::JacobiPrec, x) = y .= x ./ P.diagonal

function max_err(R)
    r = zeros(real(eltype(R)), size(R, 2))
    for j in 1:length(r)
        for i in 1:size(R, 1)
            r[j] += conj(R[i,j])*R[i,j]
        end
        r[j] = sqrt(r[j])
    end
    maximum(r)
end

@testset "Locally Optimal Block Preconditioned Conjugate Gradient" begin
    rng = Random.MersenneTwister(1234)
    @testset "Single eigenvalue" begin
        n = 10
        @testset "Small full system" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        b = rand(rng, T, n, 1)
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, largest, b; tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test norm(A*X - X*λ) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        r = lobpcg(A, largest, X; tol=10tol, log=true)
                        @test length(r.trace) == 1
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng, T, n, n)
                        B = B' + B + 20I
                        b = rand(rng, T, n, 1)
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, B, largest, b; tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*λ) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        r = lobpcg(A, B, largest, X; tol=10tol, log=true)
                        @test length(r.trace) == 1
                    end
                end
            end
        end
        @testset "Sparse Laplacian" begin
            A = laplace_matrix(Float64, 20, 2)
            rhs = randn(size(A, 2), 1)
            rmul!(rhs, inv(norm(rhs)))
            tol = IterativeSolvers.default_tolerance(Float64)
            @testset "Matrix" begin
                @testset "largest = $largest" for largest in (true, false)
                    r = lobpcg(A, largest, rhs; tol=tol, maxiter=Inf)
                    λ, X = r.λ, r.X
                    @test norm(A*X - X*λ) ≤ tol
                end
            end
        end
        @testset "Zero initial solution" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        rng_temp = Random.MersenneTwister(1234) # Issue #316 (test sensitive to the rng)
                        A = rand(rng_temp, T, n, n)
                        A = A' + A + 20I
                        b = zeros(T, n, 1)
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, largest, b; tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test norm(A*X - X*λ) ≤ tol
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        rng_temp = Random.MersenneTwister(1234) # Issue #316 (test sensitive to the rng)
                        A = rand(rng_temp, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng_temp, T, n, n)
                        B = B' + B + 20I
                        b = zeros(T, n, 1)
                        tol = IterativeSolvers.default_tolerance(T)

                        r = lobpcg(A, B, largest, b; tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*λ) ≤ tol
                    end
                end
            end
        end
        @testset "No initial solution" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        tol = IterativeSolvers.default_tolerance(T)

                        r = lobpcg(A, largest, 1; tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test norm(A*X - X*λ) ≤ tol
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng, T, n, n)
                        B = B' + B + 20I
                        tol = IterativeSolvers.default_tolerance(T)

                        r = lobpcg(A, B, largest, 1; tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*λ) ≤ tol
                    end
                end
            end
        end
        @testset "Inplace" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        tol = IterativeSolvers.default_tolerance(T)
                        b = rand(rng, T, n, 1)
                        itr = LOBPCGIterator(A, largest, b)

                        r = lobpcg!(itr; tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test norm(A*X - X*λ) ≤ tol
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng, T, n, n)
                        B = B' + B + 20I
                        b = rand(rng, T, n, 1)
                        tol = IterativeSolvers.default_tolerance(T)
                        itr = LOBPCGIterator(A, B, largest, b)

                        r = lobpcg!(itr; tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*λ) ≤ tol
                    end
                end
            end
        end
        @testset "Jacobi preconditioner" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        tol = IterativeSolvers.default_tolerance(T)
                        P = JacobiPrec(diag(A))
                        r = lobpcg(A, largest, 1; P=P, tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test norm(A*X - X*λ) ≤ tol
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        P = JacobiPrec(diag(A))
                        B = rand(rng, T, n, n)
                        B = B' + B + 20I
                        tol = IterativeSolvers.default_tolerance(T)

                        r = lobpcg(A, B, largest, 1; P=P, tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*λ) ≤ tol
                    end
                end
            end
        end
        @testset "Constraint" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, largest, 1; tol=tol, maxiter=Inf, log=false)
                        λ1, X1 = r.λ, r.X
                        r = lobpcg(A, largest, 1; C=copy(r.X), tol=tol, maxiter=Inf, log=false)
                        λ2, X2 = r.λ, r.X
                        @test norm(A*X2 - X2*λ2) ≤ tol
                        @test isapprox(real((adjoint(X1)*X2)[1,1]), 0, atol=2*n*tol)
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng, T, n, n)
                        B = B' + B + 20I
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, B, largest, 1; tol=tol, maxiter=Inf, log=false)
                        λ1, X1 = r.λ, r.X
                        r = lobpcg(A, B, largest, 1; C=copy(r.X), tol=tol, maxiter=Inf, log=false)
                        λ2, X2 = r.λ, r.X
                        @test norm(A*X2 - B*X2*λ2) ≤ tol
                        @test isapprox(real((adjoint(X1)*(B*X2))[1,1]), 0, atol=2*n*tol)
                    end
                end
            end
        end
    end
    @testset "Two eigenvalues" begin
        @testset "Small full system" begin
            n = 10
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        b = rand(rng, T, n, 2)
                        tol = IterativeSolvers.default_tolerance(T)

                        r  = lobpcg(A, largest, b; tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - X*Matrix(Diagonal(λ))) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        r = lobpcg(A, largest, X; tol=10tol, log=true)
                        @test length(r.trace) == 1
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        rng_temp = Random.MersenneTwister(123) # Issue #316 (test sensitive to the rng)
                        A = rand(rng_temp, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng_temp, T, n, n)
                        B = B' + B + 20I
                        b = rand(rng_temp, T, n, 2)
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, B, largest, b; tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*Matrix(Diagonal(λ))) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        r = lobpcg(A, B, largest, X; tol=10tol, log=true)
                        @test length(r.trace) == 1
                    end
                end
            end
        end
    end
    @testset "nev = 3, block size = $block_size" for block_size in (1, 2)
        n = 10
        @testset "Simple eigenvalue problem" begin
            @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                @testset "largest = $largest" for largest in (true, false)
                    A = rand(rng, T, n, n)
                    A = A' + A + 20I
                    tol = IterativeSolvers.default_tolerance(T)
                    X0 = rand(rng, T, n, block_size)
                    r = lobpcg(A, largest, X0, 3, tol=tol, maxiter=Inf, log=true)
                    λ, X = r.λ, r.X
                    @test max_err(A*X - X*Matrix(Diagonal(λ))) ≤ tol
                    @test all(isapprox.(adjoint(X)*X, Matrix{T}(I, 3, 3), atol=2*n*tol))
                end
            end
        end
        @testset "Generalized eigenvalue problem" begin
            @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                @testset "largest = $largest" for largest in (true, false)
                    rng_temp = Random.MersenneTwister(123) # Issue #316 (test sensitive to the rng)
                    A = rand(rng_temp, T, n, n)
                    A = A' + A + 20I
                    B = rand(rng_temp, T, n, n)
                    B = B' + B + 20I
                    tol = IterativeSolvers.default_tolerance(T)
                    X0 = rand(rng_temp, T, n, block_size)
                    r = lobpcg(A, B, largest, X0, 3, tol=tol, maxiter=Inf, log=true)
                    λ, X = r.λ, r.X
                    @test max_err(A*X - B*X*Matrix(Diagonal(λ))) ≤ tol
                    @test all(isapprox.(adjoint(X)*(B*X), Matrix{T}(I, 3, 3), atol=2*n*tol))
                end
            end
        end
        @testset "Constraint" begin
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, largest, 1; tol=tol, maxiter=Inf, log=false)
                        λ1, X1 = r.λ, r.X

                        X0 = rand(rng, T, n, block_size)
                        r = lobpcg(A, largest, X0, 3, C=copy(r.X), tol=tol, maxiter=Inf, log=true)
                        λ2, X2 = r.λ, r.X
                        @test max_err(A*X2 - X2*Matrix(Diagonal(λ2))) ≤ tol
                        @test all(isapprox.(adjoint(X2)*X2, Matrix{T}(I, 3, 3), atol=2*n*tol))
                        @test all(isapprox.(real(adjoint(X1)*X2), 0, atol=2*n*tol))
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, ComplexF32, ComplexF64)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(rng, T, n, n)
                        A = A' + A + 20I
                        B = rand(rng, T, n, n)
                        B = B' + B + 20I
                        tol = IterativeSolvers.default_tolerance(T)
                        r = lobpcg(A, B, largest, 1; tol=tol, maxiter=Inf, log=false)
                        λ1, X1 = r.λ, r.X

                        X0 = rand(rng, T, n, block_size)
                        r = lobpcg(A, B, largest, X0, 2, C=copy(r.X), tol=tol, maxiter=Inf, log=true)
                        λ2, X2 = r.λ, r.X
                        @test max_err(A*X2 - B*X2*Matrix(Diagonal(λ2))) ≤ tol
                        @test all(isapprox.(adjoint(X2)*(B*X2), Matrix{T}(I, 2, 2), atol=2*n*tol))
                        @test all(isapprox.(real(adjoint(X1)*(B*X2)), 0, atol=2*n*tol))
                    end
                end
            end
        end
    end
end

end # module
