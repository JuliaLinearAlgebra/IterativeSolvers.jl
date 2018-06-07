using IterativeSolvers
using LinearMaps
using Base.Test

include("laplace_matrix.jl")

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
    srand(1234321)
    @testset "Single eigenvalue" begin
        @testset "Small full system" begin
            n = 10
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(T, n, n)
                        A = A' * A + I
                        b = rand(T, n, 1)
                        tol = √eps(real(T))

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
                @testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(T, n, n)
                        A = A' * A + I
                        B = rand(T, n, n)
                        B = B' * B + I
                        b = rand(T, n, 1)
                        tol = √eps(real(T))

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
            scale!(rhs, inv(norm(rhs)))
            tol = 1e-5

            @testset "Matrix" begin
                @testset "largest = $largest" for largest in (true, false)
                    r = lobpcg(A, largest, rhs; tol=tol, maxiter=Inf)
                    λ, X = r.λ, r.X
                    @test norm(A*X - X*λ) ≤ tol
                end
            end
        end
    end
    @testset "Two eigenvalues" begin
        @testset "Small full system" begin
            n = 10
            @testset "Simple eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(T, n, n)
                        A = A' * A + I
                        b = rand(T, n, 2)
                        tol = √eps(real(T))

                        r  = lobpcg(A, largest, b; tol=tol, maxiter=Inf, log=false)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - X*diagm(λ)) ≤ tol
                        
                        # If you start from the exact solution, you should converge immediately
                        r = lobpcg(A, largest, X; tol=10tol, log=true)
                        @test length(r.trace) == 1
                    end
                end
            end
            @testset "Generalized eigenvalue problem" begin
                @testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
                    @testset "largest = $largest" for largest in (true, false)
                        A = rand(T, n, n)
                        A = A' * A + I
                        B = rand(T, n, n)
                        B = B' * B + I
                        b = rand(T, n, 2)
                        tol = √eps(real(T))

                        r = lobpcg(A, B, largest, b; tol=tol, maxiter=Inf, log=true)
                        λ, X = r.λ, r.X
                        @test max_err(A*X - B*X*diagm(λ)) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        r = lobpcg(A, B, largest, X; tol=10tol, log=true)
                        @test length(r.trace) == 1
                    end
                end
            end
        end
    end
end
