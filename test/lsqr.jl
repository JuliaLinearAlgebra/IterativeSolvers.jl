using IterativeSolvers

# Adapted from the BSD-licensed Matlab implementation at
#    http://www.stanford.edu/group/SOL/software/lsqr.html

# If m = n and damp = 0, this sets up a system Ax = b
# and calls lsqrSOL.m to solve it.  Otherwise, the usual
# least-squares or damped least-squares problem is solved.

# 11 Apr 1996: First version for distribution with lsqr.m.
# 07 Aug 2002: LSQR's output parameter rnorm changed to r1norm, r2norm.
# 03 May 2007: Allow A to be a matrix or a function handle.
#              Private function Aprodxxx defines matrix-vector products
#              for a specific A.
# 24 Dec 2010: A*v and A'*v use inputs (v,1) and (v,2), not (1,v) and (2,v).

#              Michael Saunders, Systems Optimization Laboratory,
#              Dept of MS&E, Stanford University.
#-----------------------------------------------------------------------
function lsqrSOLtest( m, n, damp )
    fmul = (out, b) -> Aprodxxx!(out,b,1,m,n)
    fcmul = (out, b) -> Aprodxxx!(out,b,2,m,n)
    A = MatrixCFcn{Float64}(m, n, fmul, fcmul)

    xtrue  = n : -1: 1
    b      = A*xtrue

    x = lsqr(A, b, atol = 1e-6, btol = 1e-6, conlim = 1e10, itermax = 10n)

    r    = b - A*x
    @assert norm(r) < 1e-4
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
