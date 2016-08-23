export lsmr, lsmr!

using Base.LinAlg

##############################################################################
## LSMR
##
## Minimize ||Ax-b||^2 + λ^2 ||x||^2
##
## Adapted from the BSD-licensed Matlab implementation at
## http://web.stanford.edu/group/SOL/software/lsmr/
##
## A is a StridedVecOrMat or anything that implements
## A_mul_B!(α, A, b, β, c) updates c -> α Ab + βc
## Ac_mul_B!(α, A, b, β, c) updates c -> α A'b + βc
## eltype(A)
## size(A)
## (this includes SparseMatrixCSC)
## x, v, h, hbar are AbstractVectors or anything that implements
## norm(x)
## copy!(x1, x2)
## scale!(x, α)
## axpy!(α, x1, x2)
## similar(x, T)
## length(x)
## b is an AbstractVector or anything that implements
## eltype(b)
## norm(b)
## copy!(x1, x2)
## fill!(b, α)
## scale!(b, α)
## similar(b, T)
## length(b)

##############################################################################


## Arguments:
## x is initial x0. Transformed in place to the solution.
## b equals initial b. Transformed in place
## v, h, hbar are storage arrays of length size(A, 2)
function lsmr_method!(log::ConvergenceHistory, x, A, b, v, h, hbar;
    atol::Number = 1e-6, btol::Number = 1e-6, conlim::Number = 1e8,
    maxiter::Integer = max(size(A,1), size(A,2)), λ::Number = 0)

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
    # form the first vectors u and v (satisfy  β*u = b,  α*v = A'u)
    u = A_mul_B!(-1, A, x, 1, b)
    β = norm(u)
    β > 0 && scale!(u, inv(β))
    Ac_mul_B!(1, A, u, 0, v)
    α = norm(v)
    α > 0 && scale!(v, inv(α))

    log[:atol] = atol
    log[:btol] = btol
    log[:ctol] = ctol

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
    iter = 0
    # Exit if b = 0 or A'b = 0.
    if normAr != 0
        while iter < maxiter
            nextiter!(log)
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
            push!(log, :cnorm, test3)
            push!(log, :anorm, test2)
            push!(log, :rnorm, test1)

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
    shrink!(log)
    setmvps(log, 2*iter)
    setconv(log, istop ∉ (3, 6, 7))
    x
end

## Arguments:
## x is initial x0. Transformed in place to the solution.
function lsmr!(x, A, b; maxiter::Integer = max(size(A,1)), kwargs...)
    history = ConvergenceHistory()
    reserve!(history,:anorm,maxiter)
    reserve!(history,:rnorm,maxiter)
    reserve!(history,:cnorm,maxiter)

    T = Adivtype(A, b)
    m, n = size(A, 1), size(A, 2)
    btmp = similar(b, T)
    copy!(btmp, b)
    v, h, hbar = similar(x, T), similar(x, T), similar(x, T)
    lsmr_method!(history, x, A, btmp, v, h, hbar; maxiter=maxiter, kwargs...)
    x, history
end

function lsmr(A, b; kwargs...)
    lsmr!(zerox(A, b), A, b; kwargs...)
end

for (name, symbol) in ((:Ac_mul_B!, 'T'), (:A_mul_B!, 'N'))
    @eval begin
        function Base.$name(α::Number, A::StridedVecOrMat, x::AbstractVector, β::Number, y::AbstractVector)
            BLAS.gemm!($symbol, 'N', convert(eltype(y), α), A, x, convert(eltype(y), β), y)
        end
    end
end
