export lsqr, lsqr!, master_lsqr, master_lsqr!

using Base.LinAlg # can be removed when 0.3 is no longer supported

lsqr(A, b; kwargs...) = lsqr!(zerox(A, b), A, b; kwargs...)

function lsqr!(x, A, b; kwargs...)
    lsqr_method!(x, A, b; kwargs...)
    x
end

master_lsqr(A, b; kwargs...) = master_lsqr!(zerox(A, b), A, b; kwargs...)

function master_lsqr!(x, A, b;
    atol=sqrt(eps(Adivtype(A,b))), btol=sqrt(eps(Adivtype(A,b))),
    conlim=real(one(Adivtype(A,b)))/sqrt(eps(Adivtype(A,b))),
    plot::Bool=false, maxiter::Int=max(size(A,1), size(A,2)), kwargs...
    )
    log = MethodLog(maxiter)
    add!(log,:resnorm)
    add!(log,:anorm)
    add!(log,:rnorm)
    add!(log,:cnorm)
    T = Adivtype(A, b)
    Tr = real(T)
    ctol = conlim > 0 ? convert(Tr, 1/conlim) : zero(Tr)
    conv = lsqr_method!(x, A, b;
        atol=atol, btol=btol, conlim=conlim, maxiter=maxiter, log=log, kwargs...
        )
    shrink!(log)
    plot && showplot(log)
    x, ConvergenceHistory(conv,(atol, btol, ctol),2*iters(log)+2,log)
end

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
function lsqr_method!(x, A, b;
    damp=0, atol=sqrt(eps(Adivtype(A,b))), btol=sqrt(eps(Adivtype(A,b))),
    conlim=real(one(Adivtype(A,b)))/sqrt(eps(Adivtype(A,b))),
    maxiter::Int=max(size(A,1), size(A,2)), verbose::Bool=false,
    log::MethodLog=MethodLog()
    )
    # Sanity-checking
    m = size(A,1)
    n = size(A,2)
    length(x) == n || error("x should be of length ", n)
    length(b) == m || error("b should be of length ", m)
    for i = 1:n
        isfinite(x[i]) || error("Initial guess for x must be finite")
    end

    # Initialize
    T = Adivtype(A, b)
    Tr = real(T)
    itn = istop = 0
    ctol = conlim > 0 ? convert(Tr, 1/conlim) : zero(Tr)
    Anorm = Acond = ddnorm = res2 = xnorm = xxnorm = z = sn2 = zero(Tr)
    cs2 = -one(Tr)
    dampsq = abs2(damp)
    tmpm = similar(b, T, m)
    tmpn = similar(x, T, n)

    # Set up the first vectors u and v for the bidiagonalization.
    # These satisfy  beta*u = b-A*x,  alpha*v = A'u.
    u = b - A*x
    v = copy(x)
    beta = norm(u)
    alpha = zero(Tr)
    if beta > 0
        scale!(u, inv(beta))
        Ac_mul_B!(v,A,u)
        alpha = norm(v)
    end
    if alpha > 0
        scale!(v, inv(alpha))
    end
    w = copy(v)
    wrho = similar(w)

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
        itn += 1

        # Perform the next step of the bidiagonalization to obtain the
        # next beta, u, alpha, v.  These satisfy the relations
        #      beta*u  =  A*v  - alpha*u,
        #      alpha*v  =  A'*u - beta*v.

        # Note that the following three lines are a band aid for a GEMM: X: C := αAB + βC.
        # This is already supported in A_mul_B! for sparse and distributed matrices, but not yet dense
        A_mul_B!(tmpm, A, v)
        scale!(u, -alpha)
        LinAlg.axpy!(one(eltype(tmpm)), tmpm, u)
        beta = norm(u)
        if beta > 0
            scale!(u, inv(beta))
            Anorm = sqrt(abs2(Anorm) + abs2(alpha) + abs2(beta) + dampsq)
            # Note that the following three lines are a band aid for a GEMM: X: C := αA'B + βC.
            # This is already supported in Ac_mul_B! for sparse and distributed matrices, but not yet dense
            Ac_mul_B!(tmpn, A, u)
            scale!(v, -beta)
            LinAlg.axpy!(one(eltype(tmpn)), tmpn, v)
            alpha  = norm(v)
            if alpha > 0
                scale!(v, inv(alpha))
            end
        end

        # Use a plane rotation to eliminate the damping parameter.
        # This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
        rhobar1 = sqrt(abs2(rhobar) + dampsq)
        cs1     = rhobar/rhobar1
        sn1     = damp  /rhobar1
        psi     = sn1*phibar
        phibar  = cs1*phibar

        # Use a plane rotation to eliminate the subdiagonal element (beta)
        # of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
        rho     =   sqrt(abs2(rhobar1) + abs2(beta))
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

        LinAlg.axpy!(t1, w, x)
        scale!(w, t2)
        LinAlg.axpy!(one(t2), v, w)
        copy!(wrho, w)
        scale!(wrho, inv(rho))
        ddnorm += norm(wrho)

        # Use a plane rotation on the right to eliminate the
        # super-diagonal element (theta) of the upper-bidiagonal matrix.
        # Then use the result to estimate  norm(x).
        delta   =   sn2*rho
        gambar  =  -cs2*rho
        rhs     =   phi - delta*z
        zbar    =   rhs/gambar
        xnorm   =   sqrt(xxnorm + abs2(zbar))
        gamma   =   sqrt(abs2(gambar) + abs2(theta))
        cs2     =   gambar/gamma
        sn2     =   theta /gamma
        z       =   rhs   /gamma
        xxnorm +=   abs2(z)

        # Test for convergence.
        # First, estimate the condition of the matrix  Abar,
        # and the norms of  rbar  and  Abar'rbar.
        Acond   =   Anorm*sqrt(ddnorm)
        res1    =   abs2(phibar)
        res2    =   res2 + abs2(psi)
        rnorm   =   sqrt(res1 + res2)
        Arnorm  =   alpha*abs(tau)

        next!(log)
        # 07 Aug 2002:
        # Distinguish between
        #    r1norm = ||b - Ax|| and
        #    r2norm = rnorm in current code
        #           = sqrt(r1norm^2 + damp^2*||x||^2).
        #    Estimate r1norm from
        #    r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
        # Although there is cancellation, it might be accurate enough.
        r1sq    =   abs2(rnorm) - dampsq*xxnorm
        r1norm  =   sqrt(abs(r1sq));   if r1sq < 0 r1norm = - r1norm; end
        r2norm  =   rnorm
        push!(log, :resnorm, r1norm)

        # Now use these norms to estimate certain other quantities,
        # some of which will be small near a solution.
        test1   =   rnorm /bnorm
        test2   =   Arnorm/(Anorm*rnorm)
        test3   =   inv(Acond)
        t1      =   test1/(1 + Anorm*xnorm/bnorm)
        rtol    =   btol + atol*Anorm*xnorm/bnorm

        push!(log, :cnorm, test3)
        push!(log, :anorm, test2)
        push!(log, :rnorm, test1)

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
    end
    istop > 0
end
