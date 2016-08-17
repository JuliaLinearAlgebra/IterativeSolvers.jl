export lsmr, lsmr!

using Base.LinAlg

####################
# API method calls #
####################

function lsmr(A, b; kwargs...)
    lsmr!(zerox(A, b), A, b; kwargs...)
end

function lsmr!(x, A, b;
    atol::Number = 1e-6, btol::Number = 1e-6, conlim::Number = 1e8,
    maxiter::Integer = max(size(A,1), size(A,2)), plot::Bool=false, log::Bool=false,
    kwargs...
    )
    if log
        history = ConvergenceHistory()
        history[:atol] = atol
        history[:btol] = btol
        reserve!(history,:anorm,maxiter)
        reserve!(history,:rnorm,maxiter)
        reserve!(history,:cnorm,maxiter)
    else
        history = DummyHistory()
    end
    T = Adivtype(A, b)
    m, n = size(A, 1), size(A, 2)
    btmp = similar(b, T)
    copy!(btmp, b)
    v, h, hbar = similar(x, T), similar(x, T), similar(x, T)
    lsmr_method!(x, A, btmp, v, h, hbar;
        atol=atol,btol=btol,conlim=conlim,maxiter=maxiter,log=history,kwargs...
        )
    if log
        shrink!(history)
        plot && showplot(history)
        x, history
    else
        x
    end
end

#########################
# Method Implementation #
#########################

function lsmr_method!(x, A, b, v, h, hbar;
    atol::Number = 1e-6, btol::Number = 1e-6, conlim::Number = 1e8,
    maxiter::Integer = max(size(A,1), size(A,2)), λ::Number = 0,
    log::MethodLog=DummyHistory(), verbose::Bool=false
    )
    verbose && @printf("=== lsmr ===\n%4s\t%7s\t\t%7s\t\t%7s\n","iter","anorm","cnorm","rnorm")

    # Sanity-checking
    m = size(A, 1)
    n = size(A, 2)
    length(x) == n || error("x has length $(length(x)) but should have length $n")
    length(v) == n || error("v has length $(length(v)) but should have length $n")
    length(h) == n || error("h has length $(length(h)) but should have length $n")
    length(hbar) == n || error("hbar has length $(length(hbar)) but should have length $n")
    length(b) == m || error("b has length $(length(b)) but should have length $m")

    T = Adivtype(A, b)
    Tr = real(T)
    normrs = Tr[]
    normArs = Tr[]
    conlim > 0 ? ctol = convert(Tr, inv(conlim)) : ctol = zero(Tr)
    log[:ctol] = ctol
    # form the first vectors u and v (satisfy  β*u = b,  α*v = A'u)
    u = A_mul_B!(-1, A, x, 1, b)
    β = norm(u)
    β > 0 && scale!(u, inv(β))
    Ac_mul_B!(1, A, u, 0, v)
    α = norm(v)
    α > 0 && scale!(v, inv(α))

    # Initialize variables for 1st iteration.
    ζbar = α * β
    αbar = α
    ρ = one(Tr)
    ρbar = one(Tr)
    cbar = one(Tr)
    sbar = zero(Tr)

    copy!(h, v)
    fill!(hbar, zero(Tr))

    # Initialize variables for estimation of ||r||.
    βdd = β
    βd = zero(Tr)
    ρdold = one(Tr)
    τtildeold = zero(Tr)
    θtilde  = zero(Tr)
    ζ = zero(Tr)
    d = zero(Tr)

    # Initialize variables for estimation of ||A|| and cond(A).
    normA, condA, normx = -one(Tr), -one(Tr), -one(Tr)
    normA2 = abs2(α)
    maxrbar = zero(Tr)
    minrbar = 1e100

    # Items for use in stopping rules.
    normb = β
    istop = 0
    normr = β
    normAr = α * β
    tests = Tuple{Tr, Tr, Tr}[]
    iter = 0
    # Exit if b = 0 or A'b = 0.
    if normAr != 0
        while iter < maxiter
            iter += 1
            A_mul_B!(1, A, v, -α, u)
            β = norm(u)
            if β > 0
                scale!(u, inv(β))
                Ac_mul_B!(1, A, u, -β, v)
                α = norm(v)
                α > 0 && scale!(v, inv(α))
            end

            # Construct rotation Qhat_{k,2k+1}.
            αhat = hypot(αbar, λ)
            chat = αbar / αhat
            shat = λ / αhat

            # Use a plane rotation (Q_i) to turn B_i to R_i.
            ρold = ρ
            ρ = hypot(αhat, β)
            c = αhat / ρ
            s = β / ρ
            θnew = s * α
            αbar = c * α

            # Use a plane rotation (Qbar_i) to turn R_i^T to R_i^bar.
            ρbarold = ρbar
            ζold = ζ
            θbar = sbar * ρ
            ρtemp = cbar * ρ
            ρbar = hypot(cbar * ρ, θnew)
            cbar = cbar * ρ / ρbar
            sbar = θnew / ρbar
            ζ = cbar * ζbar
            ζbar = - sbar * ζbar

            # Update h, h_hat, x.
            scale!(hbar, - θbar * ρ / (ρold * ρbarold))
            axpy!(1, h, hbar)
            axpy!(ζ / (ρ * ρbar), hbar, x)
            scale!(h, - θnew / ρ)
            axpy!(1, v, h)

            ##############################################################################
            ##
            ## Estimate of ||r||
            ##
            ##############################################################################

            # Apply rotation Qhat_{k,2k+1}.
            βacute = chat * βdd
            βcheck = - shat * βdd

            # Apply rotation Q_{k,k+1}.
            βhat = c * βacute
            βdd = - s * βacute

            # Apply rotation Qtilde_{k-1}.
            θtildeold = θtilde
            ρtildeold = hypot(ρdold, θbar)
            ctildeold = ρdold / ρtildeold
            stildeold = θbar / ρtildeold
            θtilde = stildeold * ρbar
            ρdold = ctildeold * ρbar
            βd = - stildeold * βd + ctildeold * βhat

            τtildeold = (ζold - θtildeold * τtildeold) / ρtildeold
            τd = (ζ - θtilde * τtildeold) / ρdold
            d += abs2(βcheck)
            normr = sqrt(d + abs2(βd - τd) + abs2(βdd))

            # Estimate ||A||.
            normA2 += abs2(β)
            normA  = sqrt(normA2)
            normA2 += abs2(α)

            # Estimate cond(A).
            maxrbar = max(maxrbar, ρbarold)
            if iter > 1
                minrbar = min(minrbar, ρbarold)
            end
            condA = max(maxrbar, ρtemp) / min(minrbar, ρtemp)

            ##############################################################################
            ##
            ## Test for convergence
            ##
            ##############################################################################

            # Compute norms for convergence testing.
            normAr  = abs(ζbar)
            normx = norm(x)

            # Now use these norms to estimate certain other quantities,
            # some of which will be small near a solution.
            test1 = normr / normb
            test2 = normAr / (normA * normr)
            test3 = inv(condA)
            nextiter!(log)
            push!(log, :cnorm, test3)
            push!(log, :anorm, test2)
            push!(log, :rnorm, test1)
            verbose && @printf("%3d\t%1.2e\t%1.2e\t%1.2e\n",iter,test2,test3,test1)

            t1 = test1 / (one(Tr) + normA * normx / normb)
            rtol = btol + atol * normA * normx / normb
            # The following tests guard against extremely small values of
            # atol, btol or ctol.  (The user may have set any or all of
            # the parameters atol, btol, conlim  to 0.)
            # The effect is equivalent to the normAl tests using
            # atol = eps,  btol = eps,  conlim = 1/eps.
            if iter >= maxiter istop = 7; break end
            if 1 + test3 <= 1 istop = 6; break end
            if 1 + test2 <= 1 istop = 5; break end
            if 1 + t1 <= 1 istop = 4; break end
            # Allow for tolerances set by the user.
            if test3 <= ctol istop = 3; break end
            if test2 <= atol istop = 2; break end
            if test1 <= rtol  istop = 1; break end
        end
    end
    setmvps(log, 2*iter)
    setconv(log, istop ∉ (3, 6, 7))
    verbose && @printf("\n")
end

for (name, symbol) in ((:Ac_mul_B!, 'T'), (:A_mul_B!, 'N'))
    @eval begin
        function Base.$name(α::Number, A::StridedVecOrMat, x::AbstractVector, β::Number, y::AbstractVector)
            BLAS.gemm!($symbol, 'N', convert(eltype(y), α), A, x, convert(eltype(y), β), y)
        end
    end
end

#################
# Documentation #
#################

#Initialize parameters
doc_call = """    lsmr(A, b)
"""
doc!_call = """    lsmr!(x, A, b)
"""

doc_msg = "Minimize ||Ax-b||^2 + λ^2 ||x||^2 for A*x=b.\n"
doc!_msg = "Overwrite `x`.\n\n" * doc_msg

doc_arg = ""
doc!_arg = """* `x`: initial guess, overwrite final estimation."""

doc_version = (lsmr, doc_call, doc_msg, doc_arg)
doc!_version = (lsmr!, doc!_call, doc!_msg, doc!_arg)

#Build docs
for (func, call, msg, arg) in [doc_version, doc!_version]
@doc """
$call

$msg

The method is based on the Golub-Kahan bidiagonalization process. It is
algebraically equivalent to applying MINRES to the normal equation (ATA+λ2I)x=ATb,
but has better numerical properties, especially if A is ill-conditioned.

If `log` is set to `true` is given, method will output a tuple `x, ch`. Where
`ch` is a [`ConvergenceHistory`](@ref) object. Otherwise it will only return `x`.

The `plot` attribute can only be used when `log` is set version.

**Arguments**

$arg

* `A`: linear operator.

* `b`: right hand side.

*Keywords*

* `λ::Number = 0`: lambda.

* `atol::Number = 1e-6`, `btol::Number = 1e-6`: stopping tolerances. If both are
1.0e-9 (say), the final residual norm should be accurate to about 9 digits.
(The final `x` will usually have fewer correct digits,
depending on `cond(A)` and the size of damp).

* `conlim::Number = 1e8`: stopping tolerance.  `lsmr` terminates if an estimate
of `cond(A)` exceeds conlim.  For compatible systems Ax = b,
conlim could be as large as 1.0e+12 (say).  For least-squares
problems, conlim should be less than 1.0e+8.
Maximum precision can be obtained by setting
`atol` = `btol` = `conlim` = zero, but the number of iterations
may then be excessive.

* `maxiter::Integer = min(20,length(b))`: maximum number of iterations.

* `verbose::Bool = false`: print method information.

* `log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

* `plot::Bool = false`: plot data. (Only when `log` is set)

**Output**

*`log` is `false`:*

* `x`: approximated solution.

*`log` is `true`:*

* `x`: approximated solution.

* `ch`: convergence history.

*ConvergenceHistory keys*

* `:atol` => `::Real`: atol stopping tolerance.

* `:btol` => `::Real`: btol stopping tolerance.

* `:ctol` => `::Real`: ctol stopping tolerance.

* `:anorm` => `::Real`: anorm.

* `:rnorm` => `::Real`: rnorm.

* `:cnorm` => `::Real`: cnorm.

* `:resnom` => `::Vector`: residual norm at each iteration.

""" -> func
end
