using IterativeSolvers
using Base.Test

@testset "IDR(s)" begin

n = 10
m = 6
srand(1234567)

@testset "Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
    A = rand(T, n, n) + n * eye(T, n)
    b = rand(T, n)
    tol = √eps(real(T))

    @testset "Without residual smoothing" begin
        x, history = idrs(A, b, log=true)
        @test history.isconverged
        @test norm(A * x - b) / norm(b) ≤ tol
    end

    # with smoothing
    @testset "With residual smoothing" begin
        x, history = idrs(A, b; smoothing=true, log=true)
        @test history.isconverged
        @test norm(A*x - b) / norm(b) ≤ tol
    end
end

@testset "SparseMatrixCSC{$T}" for T in (Float64, Complex128)
    A = sprand(T, n, n, 0.5) + speye(T, n)
    b = rand(T, n)
    tol = √eps(real(T))

    x, history = idrs(A, b, log=true)
    @test history.isconverged
    @test norm(A * x - b) / norm(b) ≤ tol
end
end
