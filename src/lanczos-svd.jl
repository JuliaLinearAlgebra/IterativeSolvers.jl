export svdvals_gkl

#JuliaLang/julia#12747
if VERSION <= v"0.4.0-dev+6890"
    Base.svdvals{T, Tr}(S::Base.LinAlg.SVD{T, Tr}) = (S[:S])::Vector{Tr}
end

"""
Compute the largest singular values of a matrix A using the Golub-Kahan-Lanczos
bidiagonalization method `\cite{Golub1965}` implemented in `svdvals_gkl()`.

This implementation uses full one-sided reorthogonalization as described in
`\cite{Simon2000}`.

Inputs
------

- `A`    : The matrix or matrix-like object whose truncated SVD is desired.
           `A` must support `size()`, `A*v` and `A'*v` methods.

- `nvals`: Number of singular values requested.
           Default: 6

- `v0`   : Initial guess vector.
           Default: a randomized unit vector.

Keyword arguments
-----------------

- `maxiter`: Maximum number of iterations.
             Default: the smaller dimension of A.

- `βth`    : The threshold value of β below which an invariant subspace is
             deemed to be found.
             Default: `0.1*√eps(eltype(A))`


- `σth`    : The threshold value below which a Ritz value estimate of the
             singular value is considered to be converged.
             Default: `0.1*√eps(eltype(A))`

Outputs
-------

- `converged_values`   : The requested singular values

- `convergence_history`: A `Dict{Symbol,Any}` containing the following keys:
    - `:isconverged` : Did the calculation converge?
    - `:iters`       : Number of iterations run
    - `:mvps`        : Number of matrix-vector products computed
    - `:B`           : The `Bidiagonal` matrix computed during the
                       bidiagonalization process
    - `:β`           : Norm of Lanczos iterate
    - `:ω²`          : The Frobenius norm of the difference between A and the
                       low rank approximation to A computable from the
                       currently computed singular values [Simon2000]
    - `:vals`        : Ritz values computed at each iteration
    - `:valerrs`     : Error bounds on Ritz values at each iteration

Side effects
------------

If an invariant subspace is found which smaller than the either dimension of A,
an informational message is printed and only the singular values corresponding
to this subspace are returned.

References
----------

```bibtex
@article{Golub1965,
    author = {Golub, G. and Kahan, W.},
    doi = {10.1137/0702016},
    journal = {Journal of the Society for Industrial and Applied Mathematics
        Series B Numerical Analysis},
    volume = 2,
    number = 2,
    pages = {205--224},
    title = {Calculating the Singular Values and Pseudo-Inverse of a Matrix},
    year = 1965
}


@article{Simon2000,
    author = {Simon, Horst D. and Zha, Hongyuan},
    doi = {10.1137/S1064827597327309},
    journal = {SIAM Journal on Scientific Computing},
    number = 6,
    pages = {2257--2274},
    title = {Low-Rank Matrix Approximation Using the {L}anczos Bidiagonalization
        Process with Applications},
    volume = 21,
    year = 2000
}
```
"""
function svdvals_gkl(A, nvals::Int=6, v0=randn(size(A,2));
    maxiter::Int=minimum(size(A)),
    βth::Real = 0.1*√eps(eltype(A)),
    σth::Real = 0.1*√eps(eltype(A))
    )

    m, n = size(A)
    T = eltype(A)
    Tσ= eltype(svdfact(A[1,1])[:S])
    @assert length(v0) == n
    p = v0

    αs = T[]
    βs = T[]

    α = Inf
    r = Array(T, m)
    u = Array(T, m) #Same as r
    v = Array(T, n) #Same as p

    V = Array(T, min(m, n), 0) #List of converged right vectors
    converged_values = Tσ[] #List of converged values
    converged_values_errors = T[] #List of estimated errors in converged values

    ω² = ω²₀ = vecnorm(A)^2

    convergence_history = Dict()
    convergence_history[:isconverged] = false
    convergence_history[:ω²] = [ω²]
    convergence_history[:vals] = Vector{Tσ}[]
    convergence_history[:valerrs] = Vector{Tσ}[]
    convergence_history[:mvps] = 0

    local β, k, oldσ, Δσ
    fasterror = true
    for k=1:maxiter
        #Reorthogonalize right vectors [Simon2000]
        #Here we use doubly reorthogonalized classical Gram-Schmidt
        if m ≥ n
            p -= V*(V'p)
            p -= V*(V'p)
        end

        β = norm(p)

        v[:] = p/β
        #r = A*v
        #k>1 && (r -= β*u)
        A_mul_B!(r, A, v)
        k>1 && (axpy!(-β, u, r))

        #Reorthogonalize left vectors [Simon2000]
        if m < n
            r -= V*(V'r)
            r -= V*(V'r)
        end

        α = norm(r)
        u[:] = r/α
        #p = A'u - α*v
        Ac_mul_B!(p, A, u)
        axpy!(-α, v, p)

        β = norm(p)
        push!(αs, α)
        push!(βs, β)
        convergence_history[:mvps] += 2

        #Update Simon-Zha approximation error
        ω² -= α^2
        length(βs) > 1 && (ω² -= βs[end-1]^2)
        push!(convergence_history[:ω²], ω²)

        #Compute error bars on singular values
        d = √(α*β)

        if fasterror
            σ = svdvals(Bidiagonal(αs, βs[1:end-1], false))
            Δσ = k==1 ? [d] : svdvals_error!(Δσ, oldσ, σ, d)
            fasterror = fasterror && all(isfinite(Δσ)) && all(Δσ.>√eps(Tσ))
            oldσ = σ
        end

        if !fasterror
            S = svdfact(Bidiagonal(αs, βs[1:end-1], false))
            σ = svdvals(S)
            e1= abs(S[:U][end,:])
            e2= abs(S[:Vt][:,end])
            Δσ= Tσ[d*min(e1[i], e2[i]) for i in eachindex(σ)]
        end

        #Debug svdvals_error!
        #println("Iteration $k: $(√eps())")
        #for i=1:min(k,10)
        #    println(σ[i], '\t', Δσ[i], '\t', Δσ2[i], '\t')
        #end
        #println()
        #end

        push!(convergence_history[:vals], σ)
        push!(convergence_history[:valerrs], Δσ)

        #Check number of converged values
        converged_values = Tσ[]
        for i in eachindex(σ)
            if Δσ[i] ≤ σth
                push!(converged_values, σ[i])
                push!(converged_values_errors, Δσ[i])
            end
        end

        #true: do complete reorthogonalization
        if true
            V = [V m≥n?v:u]
        end

        #If invariant subspace has been found, stop
        if β ≤ βth
            if k != minimum(size(A))
                #In exact arithmetic, Lanczos is guaranteed to find an
                #invariant subspace corresponding to the entire range of the
                #matrix. For small test matrices it is entirely possible to
                #attain this limit, so that the Krylov subspace is of rank k.
                #For other matrices this may indicate to continue with a
                #different choice of starting vector.
                info("Invariant subspace of dimension $k found")
            end
            convergence_history[:isconverged] = true
            converged_values = σ
            break
        end

        #If at least n converged values have been found, stop
        if length(converged_values) ≥ nvals
            convergence_history[:isconverged] = true
            break
        end
    end

    convergence_history[:iters] = k
    convergence_history[:B] = Bidiagonal(αs, βs[1:end-1], false)
    convergence_history[:β] = βs

    @assert issorted(converged_values, rev=true)
    converged_values, convergence_history
end

"""
Compute approximate errors in the Rayleigh-Ritz approximation to singular
values in the Golub-Kalan-Lanczos bidiagonalization method implemented in
`svdvals_gkl()`.

Warning
-------
The approximations here are only good to ~`√eps()`, since differences smaller
than that cannot be captured reliably in floating point arithmetic.

Inputs
------

- `Δθ`: Old error vector. Will get appended to.
- `τ` : Old Ritz values, sorted either in ascending or descending order.
- `θ` : Current Ritz values, sorted in the same order as `τ`
- `β` : The norm of the current residual vector

Outputs
-------

- `Δθ`: The new error vector, updated in place.

Implementation notes
--------------------

The notation used corresponds to `\cite[2/e, Corollary 7.9.2]{Parlett1980}`.

A given Ritz value ``θ`` has an associated error bound which is the product of
the norm of the residual and last entry of its associated eigenvector
`\cite[Section 13.2]{Parlett1980}`.

`\cite{Hill1992}` provides an approximation to the last component of each
eigenvector given the current Ritz values and the previous Ritz values, using
the Cauchy interlacing theorem.

This function implements the upper bound of Theorem 3 in `\cite{Hill1992}`,
adapted to the computation of error bounds on Ritz values corresponding to
singular values.

The simple modification is to square the Ritz values to get estimates of
(nonzero) eigenvalues of ``A'A`` and apply Theorem 3 of `\cite{Hill1992}` to
the Ritz values of the eigenvalues, thus deriving the final formulae used here:

```math
\Delta\theta_{j}<\beta\times
\begin{cases}
    \sqrt{\frac{\tau_{1}-\theta_{1}}{\theta_{2}-\theta_{1}}
          \frac{\tau_{1}+\theta_{1}}{\theta_{2}+\theta_{1}}}     & i=1\\
    \sqrt{\frac{\theta_{j}-\tau_{j-1}}{\theta_{j}-\theta_{j-1}}
          \frac{\tau_{j}-\theta_{j}}{\theta_{j+1}-\theta_{j}}
          \frac{\theta_{j}+\tau_{j-1}}{\theta_{j}+\theta_{j-1}}
          \frac{\tau_{j}+\theta_{j}}{\theta_{j+1}+\theta_{j}}}   & i=2,\dots,k-1\\
    \sqrt{\frac{\theta_{k}-\tau_{k-1}}{\theta_{k}-\theta_{k-1}}
          \frac{\theta_{k}+\tau_{k-1}}{\theta_{k}+\theta_{k-1}}} & i=k
\end{cases}
```

We also sneak in an `abs()` before taking the outermost square root to avoid
roundoff error when computing very small error bounds. It also has the
advantage of rendering the formulae agnostic to sort order. So long as `τ` and
`θ` are sorted the same way, the same formulae apply since reversing the sort
order merely reverses the indices in the Cauchy interlacing theorem.

References
----------

```bibtex
@book{Parlett1980,
    author = {Parlett, Beresford N},
    address = {Philadelphia, PA},
    doi = {10.1137/1.9781611971163},
    publisher = {SIAM},
    title = {The symmetric eigenvalue problem},
    year = 1980
}

@article{Hill1992,
	Author = {Hill, R O, Jr. and B N Parlett},
	Doi = {10.1137/0613019},
	Journal = {SIAM Journal on Matrix Analysis and Applications},
	Number = 1,
	Pages = {239-247},
	Title = {Refined Interlacing Properties},
	Volume = 13,
	Year = 1992
}
```
"""
function svdvals_error!{T}(
        Δθ::AbstractVector, τ::AbstractVector{T},
        θ::AbstractVector{T}, β::Real=1.0)

    n = length(θ)
    @assert (issorted(τ, rev=true) && issorted(θ, rev=true)) ||
            (issorted(τ) && issorted(θ))
    @assert length(Δθ)==length(τ)==n-1

    #Check Cauchy interlacing property
    #for i in eachindex(τ)
    #    if !(θ[i]+length(θ)^2*eps(T) ≥ τ[i] ≥ θ[i+1]-length(θ)^2*eps(T))
    #        warn("Interlacing violated $i: $(θ[i]) ≥ $(τ[i]) ≥ $(θ[i+1])")
    #        println(min(abs(τ[i] - θ[i])/eps(T), abs(θ[i+1] - τ[i])/eps(T)))
    #    end
    #end

    Δθ[1] = β*√abs((τ[1]-θ[1])*(τ[1]+θ[1])/
               ((θ[2]-θ[1])*(θ[2]+θ[1])))

    for j=2:n-1
        Δθ[j] = β*√abs(((θ[j]-τ[j-1])/(θ[j]-θ[j-1])) *
                    ((θ[j]+τ[j-1])/(θ[j]+θ[j-1])) *
                    ((τ[j]-θ[j])/(θ[j+1]-θ[j])) *
                    ((τ[j]+θ[j])/(θ[j+1]+θ[j])))
    end
    if n>1
        push!(Δθ, β*√abs((θ[n]-τ[n-1])*(θ[n]+τ[n-1])/
                        ((θ[n]-θ[n-1])*(θ[n]+θ[n-1]))))
    end
    Δθ
end

