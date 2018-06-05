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

                        λ, x  = lobpcg(A, largest, b; tol=tol, maxiter=Inf, log=false)
                        @test norm(A*x - x*λ) ≤ tol
                        
                        # If you start from the exact solution, you should converge immediately
                        λ, x, ch = lobpcg(A, largest, x; tol=10tol, log=true)
                        @test length(ch) ≤ 1
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

                        λ, x, ch = lobpcg(A, B, largest, b; tol=tol, maxiter=Inf, log=true)
                        @show max_err(A*x - B*x*diagm(λ)), tol
                        @test max_err(A*x - B*x*diagm(λ)) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        λ, x, ch = lobpcg(A, B, largest, x; tol=10tol, log=true)
                        @show length(ch)
                        @test length(ch) ≤ 1
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
                    λ, xLOBPCG = lobpcg(A, largest, rhs; tol=tol, maxiter=Inf)
                    @test norm(A * xLOBPCG - xLOBPCG * λ) ≤ tol
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

                        λ, x  = lobpcg(A, largest, b; tol=tol, maxiter=Inf, log=false)
                        @test max_err(A*x - x*diagm(λ)) ≤ tol
                        
                        # If you start from the exact solution, you should converge immediately
                        λ, x, ch = lobpcg(A, largest, x; tol=10tol, log=true)
                        @test length(ch) ≤ 1
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

                        λ, x, ch = lobpcg(A, B, largest, b; tol=tol, maxiter=Inf, log=true)
                        @show max_err(A*x - B*x*diagm(λ)), tol
                        @test max_err(A*x - B*x*diagm(λ)) ≤ tol

                        # If you start from the exact solution, you should converge immediately
                        λ, x, ch = lobpcg(A, B, largest, x; tol=10tol, log=true)
                        @show length(ch)
                        @test length(ch) ≤ 1
                    end
                end
            end
        end
    end
end
