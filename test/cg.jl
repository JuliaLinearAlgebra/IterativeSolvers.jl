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
    x,ch = cg(A,rhs;tol=tol, maxiter=2*N)

    @fact norm(A*x - rhs) --> less_than(cond(A)*√tol)
    @fact ch.isconverged --> true

    # If you start from the exact solution, you should converge immediately
    x2,ch2 = cg!(A\rhs, A, rhs; tol=tol*10)
    @fact length(ch2.residuals) --> less_than_or_equal(1)

    # Test with cholfact should converge immediately
    F = cholfact(A)
    x2,ch2 = cg(A, rhs, F)
    @fact length(ch2.residuals) --> less_than_or_equal(2)
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
        xCG, = cg(A,rhs;tol=tol,maxiter=100)
        xJAC, = cg(A,rhs,JAC;tol=tol,maxiter=100)
        xSGS, = cg(A,rhs,SGS;tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)
    end

    Af = FuncMatrix(A)
    context("function") do
        xCG, = cg(Af,rhs;tol=tol,maxiter=100)
        xJAC, = cg(Af,rhs,JAC;tol=tol,maxiter=100)
        xSGS, = cg(Af,rhs,SGS;tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)
    end

    context("function with specified starting guess") do
        tol = 1e-4
        x0 = randn(size(rhs))
        xCG, hCG = cg!(copy(x0),Af,rhs,identity;tol=tol,maxiter=100)
        xJAC, hJAC = cg!(copy(x0),Af,rhs,JAC;tol=tol,maxiter=100)
        xSGS, hSGS = cg!(copy(x0),Af,rhs,SGS;tol=tol,maxiter=100)
        @fact norm(A*xCG - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xSGS - rhs) --> less_than_or_equal(tol)
        @fact norm(A*xJAC - rhs) --> less_than_or_equal(tol)

        iterCG = length(hCG.residuals)
        iterJAC = length(hJAC.residuals)
        iterSGS = length(hSGS.residuals)
        @fact iterJAC --> iterCG
        @fact iterSGS --> less_than_or_equal(iterJAC) "Preconditioner increased the number of iterations"
    end
end

end
