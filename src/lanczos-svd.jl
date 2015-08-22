export svdvals_gkl

#JuliaLang/julia#12747
if VERSION <= v"0.4.0-dev+6890"
    Base.svdvals{T, Tr}(S::Base.LinAlg.SVD{T, Tr}) = (S[:S])::Vector{Tr}
end

"""
Compute the largest singular values of a matrix A using the Golub-Kahan-Lanczos
bidiagonalization method [Golub1965].

This implementation uses full one-sided reorthogonalization as described in
[Simon2000].

Inputs

    A: The matrix or matrix-like object whose SVD is desired.
       A must support size(), A*v and A'*v methods.

    nvals: Number of singular values requested. (Default: 6)

    v0: Initial guess vector. Default: a randomized unit vector.


Keyword arguments

    maxiter: Maximum number of iterations. Default: the smaller dimension of A.

    βth    : The threshold value of β below which an invariant subspace is
             deemed to be found. Default: 0.1*√eps(eltype(A))

    σth    : The threshold value below which a Ritz value estimate of the
             singular value is considered to be converged. Default:
             0.1*√eps(eltype(A))


Outputs

    converged_values: The requested singular values
    B: The bidiagonal matrix produced during Step 1 of the process.

Side effects

    Prints a convergence summary consisting of:

    - The number of iterations

    - The final approximation error ω², which is the Frobenius norm of the
    difference between A and the low rank approximation to A computable from
    the computed singular values [Simon2000]

    If an invariant subspace is found which smaller than the either dimension
    of A, an informational message is printed and only the singular values
    corresponding to this subspace are returned.

References

@article{Golub1965,
    author = {Golub, G. and Kahan, W.},
    doi = {10.1137/0702016},
    journal = {Journal of the Society for Industrial and Applied Mathematics Series B Numerical Analysis},
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
    title = {Low-Rank Matrix Approximation Using the {L}anczos Bidiagonalization Process with Applications},
    volume = 21,
    year = 2000
}
"""
function svdvals_gkl(A, nvals::Int=6, v0=randn(size(A,2)), maxiter::Int=minimum(size(A)),
        βth::Real = 0.1*√eps(eltype(A)), σth::Real = 0.1*√eps(eltype(A)))

    m, n = size(A)
    T = eltype(A)
    Tσ= eltype(svdfact(A[1,1])[:S])
    p = v0

    αs = T[]
    βs = T[]

    β = norm(p)
    α = Inf
    u = T[]

    converged_vectors = Vector{T}[] #List of converged right vectors
    converged_values = Tσ[] #List of converged values
    converged_values_errors = T[] #List of estimated errors in converged values

    ω² = ω²₀ = vecnorm(A)^2

    k = 0
    for k=1:maxiter
        #Reorthogonalize right vectors [Simon2000]
        if m >= n
            for w in converged_vectors
                p -= (p⋅w)*w
            end
        end

        v = scale!(p, inv(β))
        r = A*p
        k>1 && (r -= β*u)

        #Reorthogonalize left vectors [Simon2000]
        if m < n
            for w in converged_vectors
                r -= (r⋅w)*w
            end
        end

        α = norm(r)
        u = scale!(r, inv(α))
        p = A'u - α*v
        β = norm(p)
        push!(αs, α)
        push!(βs, β)

        #Update Simon-Zha approximation error
        ω² -= α^2
        length(βs) > 1 && (ω² -= βs[end-1]^2)

        #Compute error bars on singular values
        S = svdfact(Bidiagonal(αs, βs[1:end-1], false))
        d = √(α*β)
        e1 = abs(S[:U][end, :])
        e2 = abs(S[:Vt][:, end])
        σ = svdvals(S)
        Δσ = [min(d*e1[i], d*e2[i]) for i in eachindex(e1)]

        #Check number of converged values
        converged_values = Tσ[]
        for i=1:length(αs)
            if Δσ[i] ≤ βth
                push!(converged_values, σ[i])
                push!(converged_values_errors, Δσ[i])
            end
        end

        #true: do complete reorthogonalization
        if true
            push!(converged_vectors, m>=n ? v : u)
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
            converged_values = σ
            break
        end

        #If at least n converged values have been found, stop
        length(converged_values) ≥ nvals && break
    end

    info("Convergence summary")
    info("Number of iterations: $k")
    info("Final approximation error: ω² = $(ω²) ($(100ω²/ω²₀)%)")
    info("Final Lanczos β = $β")
    sort!(converged_values, rev=true), Bidiagonal(αs, βs[1:end-1], false)
end
