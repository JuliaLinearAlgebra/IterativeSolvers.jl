export lsqr, lsqr!

####################
# API method calls #
####################

lsqr(A, b; kwargs...) = lsqr!(zerox(A, b), A, b; kwargs...)

function lsqr!(x, A, b;
    plot::Bool=false, maxiter::Int=max(size(A,1), size(A,2)), log::Bool=false,
    kwargs...
    )
    T = Adivtype(A, b)
    z = zero(T)

    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    reserve!(history,[:resnorm,:anorm,:rnorm,:cnorm],maxiter)
    lsqr_method!(history, x, A, b; maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

# Michael Saunders, Systems Optimization Laboratory,
# Dept of MS&E, Stanford University.
#
# Adapted for Julia by Timothy E. Holy with the following changes:
#    - Allow an initial guess for x
#    - Eliminate printing
#----------------------------------------------------------------------
function lsqr_method!(log::ConvergenceHistory, x, A, b;
    damp=0, atol=sqrt(eps(real(Adivtype(A,b)))), btol=sqrt(eps(real(Adivtype(A,b)))),
    conlim=real(one(Adivtype(A,b)))/sqrt(eps(real(Adivtype(A,b)))),
    maxiter::Int=max(size(A,1), size(A,2)), verbose::Bool=false,
    )
    verbose && @printf("=== lsqr ===\n%4s\t%7s\t\t%7s\t\t%7s\t\t%7s\n","iter","resnorm","anorm","cnorm","rnorm")
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

    log[:atol] = atol
    log[:btol] = btol
    log[:ctol] = ctol

    # Set up the first vectors u and v for the bidiagonalization.
    # These satisfy  beta*u = b-A*x,  alpha*v = A'u.
    u = b - A*x
    v = copy(x)
    beta = norm(u)
    alpha = zero(Tr)
    if beta > 0
        log.mtvps=1
        @blas! u *= inv(beta)
        Ac_mul_B!(v,A,u)
        alpha = norm(v)
    end
    if alpha > 0
        @blas! v *= inv(alpha)
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
    while (itn < maxiter) & !log.isconverged
        nextiter!(log,mvps=1)
        itn += 1

        # Perform the next step of the bidiagonalization to obtain the
        # next beta, u, alpha, v.  These satisfy the relations
        #      beta*u  =  A*v  - alpha*u,
        #      alpha*v  =  A'*u - beta*v.

        # Note that the following three lines are a band aid for a GEMM: X: C := αAB + βC.
        # This is already supported in A_mul_B! for sparse and distributed matrices, but not yet dense
        A_mul_B!(tmpm, A, v)
        @blas! u *= -alpha
        @blas! u += one(eltype(tmpm))*tmpm
        beta = norm(u)
        if beta > 0
            log.mtvps+=1
            @blas! u *= inv(beta)
            Anorm = sqrt(abs2(Anorm) + abs2(alpha) + abs2(beta) + dampsq)
            # Note that the following three lines are a band aid for a GEMM: X: C := αA'B + βC.
            # This is already supported in Ac_mul_B! for sparse and distributed matrices, but not yet dense
            Ac_mul_B!(tmpn, A, u)
            @blas! v *= -beta
            @blas! v += one(eltype(tmpn))*tmpn
            alpha  = norm(v)
            if alpha > 0
                @blas! v *= inv(alpha)
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

        @blas! x += t1*w
        @blas! w *= t2
        @blas! w += one(t2)*v
        @blas! wrho = w
        @blas! wrho *= inv(rho)
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
        verbose && @printf("%3d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\n",itn,r1norm,test2,test3,test1)

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

        setconv(log, istop > 0)
    end
    verbose && @printf("\n")
    x
end

#################
# Documentation #
#################

let
#Initialize parameters
doc_call = """    lsqr(A, b)
"""
doc!_call = """    lsqr!(x, A, b)
"""

doc_msg = """LSQR solves Ax = b or min ||b - Ax||^2 if damp = 0,
or   min ||(b) - (  A   )x||   otherwise.
         ||(0)   (damp*I) ||^2.
"""
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """`x`: initial guess, overwrite final estimation."""

doc_version = (lsqr, doc_call, doc_msg, doc_arg)
doc!_version = (lsqr!, doc!_call, doc!_msg, doc!_arg)

i=0
docstring = Vector(2)

#Build docs
for (func, call, msg, arg) in [doc_version, doc!_version]
i+=1
docstring[i] = """
$call

$msg

The method is based on the Golub-Kahan bidiagonalization process.
It is algebraically equivalent to applying CG to the normal equation (ATA+λ2I)x=ATb,
but has better numerical properties, especially if A is ill-conditioned.

If `log` is set to `true` is given, method will output a tuple `x, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `x`.

The `plot` attribute can only be used when `log` is set version.

# Arguments

$arg

`A`: linear operator.

`b`: right hand side.

## Keywords

`damp::Number = 0`: damping parameter.

`atol::Number = 1e-6`, `btol::Number = 1e-6`: stopping tolerances. If both are
1.0e-9 (say), the final residual norm should be accurate to about 9 digits.
(The final `x` will usually have fewer correct digits,
depending on `cond(A)` and the size of damp).

`conlim::Number = 1e8`: stopping tolerance.  `lsmr` terminates if an estimate
of `cond(A)` exceeds conlim.  For compatible systems Ax = b,
conlim could be as large as 1.0e+12 (say).  For least-squares
problems, conlim should be less than 1.0e+8.
Maximum precision can be obtained by setting
`atol` = `btol` = `conlim` = zero, but the number of iterations
may then be excessive.

`maxiter::Integer = min(20,length(b))`: maximum number of iterations.

`verbose::Bool = false`: print method information.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

`plot::Bool = false`: plot data. (Only when `log` is set)

# Output

**if `log` is `false`**

`x`: approximated solution.

**if `log` is `true`**

`x`: approximated solution.

`ch`: convergence history.

**ConvergenceHistory keys**

`:atol` => `::Real`: atol stopping tolerance.

`:btol` => `::Real`: btol stopping tolerance.

`:ctol` => `::Real`: ctol stopping tolerance.

`:anorm` => `::Real`: anorm.

`:rnorm` => `::Real`: rnorm.

`:cnorm` => `::Real`: cnorm.

`:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> lsqr
@doc docstring[2] -> lsqr!
end
