export lsqr, lsqr!

# Adapted from the BSD-licensed Matlab implementation at
#    http://www.stanford.edu/group/SOL/software/lsqr.html
#
# LSQR solves Ax = b or min ||b - Ax||_2 if damp = 0,
# or   min ||(b) - (  A   )x||   otherwise.
#          ||(0)   (damp*I) ||_2

#-----------------------------------------------------------------------
# LSQR uses an iterative (conjugate-gradient-like) method.
# For further information, see
# 1. C. C. Paige and M. A. Saunders (1982a).
#    LSQR: An algorithm for sparse linear equations and sparse least squares,
#    ACM TOMS 8(1), 43-71.
# 2. C. C. Paige and M. A. Saunders (1982b).
#    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,
#    ACM TOMS 8(2), 195-209.
# 3. M. A. Saunders (1995).  Solution of sparse rectangular systems using
#    LSQR and CRAIG, BIT 35, 588-604.
#
# Input parameters:
# atol, btol  are stopping tolerances.  If both are 1.0e-9 (say),
#             the final residual norm should be accurate to about 9 digits.
#             (The final x will usually have fewer correct digits,
#             depending on cond(A) and the size of damp.)
# conlim      is also a stopping tolerance.  lsqr terminates if an estimate
#             of cond(A) exceeds conlim.  For compatible systems Ax = b,
#             conlim could be as large as 1.0e+12 (say).  For least-squares
#             problems, conlim should be less than 1.0e+8.
#             Maximum precision can be obtained by setting
#             atol = btol = conlim = zero, but the number of iterations
#             may then be excessive.
# maxiter      is an explicit limit on iterations (for safety).
#
#              Michael Saunders, Systems Optimization Laboratory,
#              Dept of MS&E, Stanford University.
#
# Adapted for Julia by Timothy E. Holy with the following changes:
#    - Allow an initial guess for x
#    - Eliminate printing
#-----------------------------------------------------------------------
function lsqr!(x, ch::ConvergenceHistory, A, b; damp=0, atol=sqrt(eps(Abtype(A,b))), btol=sqrt(eps(Abtype(A,b))), conlim=one(Abtype(A,b))/sqrt(eps(Abtype(A,b))), maxiter::Int=max(size(A,1), size(A,2)))
    # Sanity-checking
    m = size(A,1)
    n = size(A,2)
    length(x) == n || error("x should be of length ", n)
    length(b) == m || error("b should be of length ", m)
    for i = 1:n
        isfinite(x[i]) || error("Initial guess for x must be finite")
    end

    # Initialize
    empty!(ch)
    ch.threshold = (atol, btol, conlim)
    T = Abtype(A, b)
    itn = istop = 0
    ctol = conlim > 0 ? convert(T,1/conlim) : zero(T)
    Anorm = Acond = ddnorm = res2 = xnorm = xxnorm = z = sn2 = zero(T)
    cs2 = -one(T)
    dampsq = damp*damp
    tmpm = Array(T, m)
    tmpn = Array(T, n)

    # Set up the first vectors u and v for the bidiagonalization.
    # These satisfy  beta*u = b-A*x,  alpha*v = A'u.
    u = b-A*x
    v = zeros(T, n)
    beta = norm(u)
    alpha = zero(T)
    if beta > 0
        scale!(u, one(T)/beta)
        Ac_mul_B!(v,A,u)
        alpha = norm(v)
    end
    if alpha > zero(T)
        scale!(v, one(T)/alpha)
    end
    w = copy(v)
    ch.mvps += 2

    Arnorm = alpha*beta
    if Arnorm == 0
        return
    end

    rhobar = alpha
    phibar = bnorm = rnorm = r1norm = r2norm = beta

    #------------------------------------------------------------------
    #     Main iteration loop.
    #------------------------------------------------------------------
    while itn < maxiter
        itn = itn + 1

        # Perform the next step of the bidiagonalization to obtain the
        # next beta, u, alpha, v.  These satisfy the relations
        #      beta*u  =  A*v  - alpha*u,
        #      alpha*v  =  A'*u - beta*v.
        A_mul_B!(tmpm, A, v)
        for i = 1:m
            u[i] = tmpm[i] - alpha*u[i]
        end
        beta = norm(u)
        if beta > 0
            scale!(u, one(T)/beta)
            Anorm = sqrt(Anorm*Anorm + alpha*alpha + beta*beta + dampsq)
            Ac_mul_B!(tmpn, A, u)
            for i = 1:n
                v[i] = tmpn[i] - beta*v[i]
            end
            alpha  = norm(v)
            if alpha > 0
                for i = 1:n v[i] /= alpha; end
            end
        end
        ch.mvps += 2

        # Use a plane rotation to eliminate the damping parameter.
        # This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
        rhobar1 = sqrt(rhobar*rhobar + dampsq)
        cs1     = rhobar/rhobar1
        sn1     = damp  /rhobar1
        psi     = sn1*phibar
        phibar  = cs1*phibar

        # Use a plane rotation to eliminate the subdiagonal element (beta)
        # of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
        rho     =   sqrt(rhobar1*rhobar1 + beta*beta)
        cs      =   rhobar1/rho
        sn      =   beta   /rho
        theta   =   sn*alpha
        rhobar  = - cs*alpha
        phi     =   cs*phibar
        phibar  =   sn*phibar
        tau     =   sn*phi

        # Update x and w
        t1      =   phi  /rho
        t2      = - theta/rho
        for i = 1:n
            wi = w[i]
            x[i] += t1*wi
            w[i] = v[i] + t2*wi
            wirho = wi/rho
            ddnorm += wirho*wirho
        end

        # Use a plane rotation on the right to eliminate the
        # super-diagonal element (theta) of the upper-bidiagonal matrix.
        # Then use the result to estimate  norm(x).
        delta   =   sn2*rho
        gambar  = - cs2*rho
        rhs     =   phi - delta*z
        zbar    =   rhs/gambar
        xnorm   =   sqrt(xxnorm + zbar^2)
        gamma   =   sqrt(gambar*gambar + theta*theta)
        cs2     =   gambar/gamma
        sn2     =   theta /gamma
        z       =   rhs   /gamma
        xxnorm  =   xxnorm + z^2

        # Test for convergence.
        # First, estimate the condition of the matrix  Abar,
        # and the norms of  rbar  and  Abar'rbar.
        Acond   =   Anorm*sqrt(ddnorm)
        res1    =   phibar^2
        res2    =   res2 + psi^2
        rnorm   =   sqrt(res1 + res2)
        Arnorm  =   alpha*abs(tau)

        # 07 Aug 2002:
        # Distinguish between
        #    r1norm = ||b - Ax|| and
        #    r2norm = rnorm in current code
        #           = sqrt(r1norm^2 + damp^2*||x||^2).
        #    Estimate r1norm from
        #    r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
        # Although there is cancellation, it might be accurate enough.
        r1sq    =   rnorm^2 - dampsq*xxnorm
        r1norm  =   sqrt(abs(r1sq));   if r1sq < 0 r1norm = - r1norm; end
        r2norm  =   rnorm
        push!(ch, r1norm)

        # Now use these norms to estimate certain other quantities,
        # some of which will be small near a solution.
        test1   =   rnorm /bnorm
        test2   =   Arnorm/(Anorm*rnorm)
        test3   =   one(T)/Acond
        t1      =   test1/(one(T) + Anorm*xnorm/bnorm)
        rtol    =   btol + atol*Anorm*xnorm/bnorm

        # The following tests guard against extremely small values of
        # atol, btol  or  ctol.  (The user may have set any or all of
        # the parameters  atol, btol, conlim  to 0.)
        # The effect is equivalent to the normal tests using
        # atol = eps,  btol = eps,  conlim = 1/eps.
        if itn >= maxiter  istop = 7; end
        if 1 + test3  <= 1 istop = 6; end
        if 1 + test2  <= 1 istop = 5; end
        if 1 + t1     <= 1 istop = 4; end

        # Allow for tolerances set by the user
        if  test3 <= ctol  istop = 3; end
        if  test2 <= atol  istop = 2; end
        if  test1 <= rtol  istop = 1; end

        if istop > 0 break end
    end
    x
end

function lsqr(A, b; kwargs...)
    T = Abtype(A, b)
    z = zero(T)
    x = zeros(T, size(A, 2))
    ch = ConvergenceHistory(false, (z,z,z), 0, T[])
    lsqr!(x, ch, A, b; kwargs...)
    x, ch
end
