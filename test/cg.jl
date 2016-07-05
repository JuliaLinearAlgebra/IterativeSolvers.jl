using IterativeSolvers
using FactCheck

include("getDivGrad.jl")

facts("cg") do

context("Small full system") do
    N=10
    A = randn(N,N)
    A = A'*A
    rhs = randn(N)
    tol = 1e-12
    x,ch = master_cg(A,rhs;tol=tol, maxiter=2*N)

    @fact norm(A*x - rhs) --> less_than(cond(A)*√tol)
    @fact ch.isconverged --> true

    # If you start from the exact solution, you should converge immediately
    x2,ch2 = master_cg!(A\rhs, A, rhs; tol=tol*10)
    @fact iters(ch2) --> less_than_or_equal(1)

    # Test with cholfact should converge immediately
    F = cholfact(A)
    x2,ch2 = master_cg(A, rhs; pl=F)
    @fact iters(ch2) --> less_than_or_equal(2)
end

context("Sparse Laplacian") do
    A = getDivGrad(32,32,32)
    L = tril(A)
    D = diag(A)
    U = triu(A)
    JAC(x) = D.\x
    SGS(x) = L\(D.*(U\x))

    rhs = randn(size(A,2))
    rhs/= norm(rhs)
    tol = 1e-5

    context("matrix") do
        xCG, = master_cg(A,rhs;tol=tol,maxiter=100)
        xJAC, = master_cg(A,rhs;pl=JAC,tol=tol,maxiter=100)
        xSGS, = master_cg(A,rhs;pl=SGS,tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)
    end

    Af = MatrixFcn(A)
    context("function") do
        xCG, = master_cg(Af,rhs;tol=tol,maxiter=100)
        xJAC, = master_cg(Af,rhs;pl=JAC,tol=tol,maxiter=100)
        xSGS, = master_cg(Af,rhs;pl=SGS,tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)
    end

    context("function with specified starting guess") do
        tol = 1e-4
        x0 = randn(size(rhs))
        xCG, hCG = master_cg!(copy(x0),Af,rhs;pl=identity,tol=tol,maxiter=100)
        xJAC, hJAC = master_cg!(copy(x0),Af,rhs;pl=JAC,tol=tol,maxiter=100)
        xSGS, hSGS = master_cg!(copy(x0),Af,rhs;pl=SGS,tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)

        iterCG = iters(hCG)
        iterJAC = iters(hJAC)
        iterSGS = iters(hSGS)
        @fact iterJAC --> iterCG
        @fact iterSGS --> less_than_or_equal(iterJAC) "Preconditioner increased the number of iterations"
    end
end

end
