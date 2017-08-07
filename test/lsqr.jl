using IterativeSolvers
using Base.Test
using LinearMaps

srand(1234321)

@testset "LSQR" begin

@testset "Small dense matrix" for T = (Float32, Float64)
    A = rand(T, 10, 5)
    b = rand(T, 10)
    x, history = lsqr(A, b, log = true)
    @test norm(x - A\b) ≤ √eps(T)
    @test history.isconverged
    @test last(history[:resnorm]) ≈ norm(b - A * x) atol=√eps(T)
end

function sol_matrix(m, n)
    mn = min(m, n)
    spdiagm((1.0 : mn - 1, 1.0 : mn), (-1, 0), m, n)
end

@testset "SOL test" for (m, n) = ((10, 10), (20, 10), (20, 10))
    # Test adapted from the BSD-licensed Matlab implementation at
    #    http://www.stanford.edu/group/SOL/software/lsqr.html
    #              Michael Saunders, Systems Optimization Laboratory,
    #              Dept of MS&E, Stanford University.
    #-----------------------------------------------------------------------
    A = LinearMap(sol_matrix(m, n))
    x = float.(n : -1 : 1)
    b = float(A*x)
    x_lsqr = lsqr(A, b, atol=1e-6, btol=1e-6, conlim=1e10, maxiter=10n)
    @test norm(b - A * x_lsqr) ≤ 1e-4
end
end
