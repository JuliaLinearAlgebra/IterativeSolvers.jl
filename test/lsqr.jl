using IterativeSolvers
using FactCheck
using Base.Test

facts("lsqr") do

context("Small dense matrix") do
    A = rand(10, 5)
    b = rand(10)
    x = lsqr(A, b)
    @fact norm(x - A\b) --> less_than(√eps())
end

context("SOL test") do
    # Test adapted from the BSD-licensed Matlab implementation at
    #    http://www.stanford.edu/group/SOL/software/lsqr.html
    #              Michael Saunders, Systems Optimization Laboratory,
    #              Dept of MS&E, Stanford University.
    #-----------------------------------------------------------------------
    function lsqrSOLtest( m, n, damp )
        fmul  = (out, b) -> Aprodxxx!(out,b,1,m,n)
        fcmul = (out, b) -> Aprodxxx!(out,b,2,m,n)
        A = MatrixCFcn{Int}(m, n, fmul, fcmul)
        xtrue = n:-1:1
        b = float(A*xtrue)
        x = lsqr(A, b, atol = 1e-6, btol = 1e-6, conlim = 1e10, maxiter = 10n)
        r = b - A*x
        @fact norm(r) --> less_than_or_equal(1e-4)
        x
    end

    function Aprodxxx!(y, x, mode, m, n )
        # if mode = 1, computes y = A*x
        # if mode = 2, computes y = A'*x
        # for some matrix  A.
        #
        # This is a simple example for testing  LSQR.
        # It uses the leading m*n submatrix from
        # A = [ 1
        #       1 2
        #         2 3
        #           3 4
        #             ...
        #               n ]
        # suitably padded by zeros.
        #
        # 11 Apr 1996: First version for distribution with lsqr.m.
        #              Michael Saunders, Dept of EESOR, Stanford University.

        if mode == 1
            y[1] = x[1]
            for i = 2:n
                y[i] = i*x[i] + (i-1)*x[i-1]
            end
            for i = n+1:m
                y[i] = 0
            end
        else
            mn = min(m, n)
            for i = 1:mn-1
                y[i] = i*(x[i]+x[i+1])
            end
            y[mn] = mn*x[mn]
            for i = m+1:n
                y[i] = 0
            end
        end
        y
    end

    x = lsqrSOLtest( 10,10, 0   )
    x = lsqrSOLtest( 20,10, 0   )
    x = lsqrSOLtest( 20,10, 0.1 )
end

context("Issue 64") do
    srand(224)
    n=10
    A=sprand(n,n,.5)
    b=rand(n)

    x, ch = lsqr(A,b,maxiter=100, log=true)
    resnorm = norm(A*x - b)
    @fact resnorm --> less_than(√eps())
    @fact ch.isconverged --> true
    @fact last(ch[:resnorm]) --> roughly(resnorm, atol=√eps())
end

end
