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
    x,ch = cg(Master, A,rhs;tol=tol, maxiter=2*N)

    @fact norm(A*x - rhs) --> less_than(cond(A)*√tol)
    @fact ch.isconverged --> true

    # If you start from the exact solution, you should converge immediately
    x2,ch2 = cg!(Master, A\rhs, A, rhs; tol=tol*10)
    @fact niters(ch2) --> less_than_or_equal(1)

    # Test with cholfact should converge immediately
    F = cholfact(A)
    x2,ch2 = cg(Master, A, rhs; pl=F)
    @fact niters(ch2) --> less_than_or_equal(2)
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
        xCG, = cg(Master, A,rhs;tol=tol,maxiter=100)
        xJAC, = cg(Master, A,rhs;pl=JAC,tol=tol,maxiter=100)
        xSGS, = cg(Master, A,rhs;pl=SGS,tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)
    end

    Af = MatrixFcn(A)
    context("function") do
        xCG, = cg(Master, Af,rhs;tol=tol,maxiter=100)
        xJAC, = cg(Master, Af,rhs;pl=JAC,tol=tol,maxiter=100)
        xSGS, = cg(Master, Af,rhs;pl=SGS,tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)
    end

    context("function with specified starting guess") do
        tol = 1e-4
        x0 = randn(size(rhs))
        xCG, hCG = cg!(Master, copy(x0),Af,rhs;pl=identity,tol=tol,maxiter=100)
        xJAC, hJAC = cg!(Master, copy(x0),Af,rhs;pl=JAC,tol=tol,maxiter=100)
        xSGS, hSGS = cg!(Master, copy(x0),Af,rhs;pl=SGS,tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)

        iterCG = niters(hCG)
        iterJAC = niters(hJAC)
        iterSGS = niters(hSGS)
        @fact iterJAC --> iterCG
        @fact iterSGS --> less_than_or_equal(iterJAC) "Preconditioner increased the number of iterations"
    end
end

end
