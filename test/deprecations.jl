module TestDeprecations

using IterativeSolvers
using Test

@testset "Test depwarn for tol" begin
    A = rand(2, 2)
    b = rand(2)

    @test_deprecated cg(A, b, tol=1.0, maxiter=1)
    @test_deprecated bicgstabl(A, b, 1, tol=1.0, max_mv_products=1)
    @test_deprecated chebyshev(A, b, 0.0, 1.0, tol=1.0, maxiter=1)
    @test_deprecated gmres(A, b, tol=1.0, maxiter=1)
    @test_deprecated idrs(A, b, tol=1.0, maxiter=1)
    @test_deprecated minres(A, b, tol=1.0, maxiter=1)
    @test_deprecated qmr(A, b, tol=1.0, maxiter=1)
end

end # module
