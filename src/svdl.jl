import Base.LinAlg: axpy!

export svdl

"""
    BrokenArrowBidiagonal

Matrix of the form
 d1          a_1
    d2       a_2
      ...    ...
         d_l a_l
             d_l+1 e_l+1
                   ...  e_k-1
                        d_k

# Implements

* Base: size, getindex, full, svdfact

"""
type BrokenArrowBidiagonal{T} <: AbstractMatrix{T}
    dv::Vector{T}
    av::Vector{T}
    ev::Vector{T}
end

Base.size(B::BrokenArrowBidiagonal) = (n=length(B.dv); (n, n))
function Base.size(B::BrokenArrowBidiagonal, n::Int)
    if n==1 || n==2
        return length(B.dv)
    else
        throw(ArgumentError("invalid dimension $n"))
    end
end

function Base.getindex{T}(B::BrokenArrowBidiagonal{T}, i::Int, j::Int)
    n = size(B, 1)
    k = length(B.av)
    if !(1 ≤ i ≤ n || 1 ≤ j ≤ n)
        throw(BoundsError())
    end

    if i == j
        return B.dv[i]
    elseif i ≤ k && j == k+1
        return B.av[i]
    elseif i > k && j == i+1
        return B.ev[i-k]
    else
        return zero(T)
    end
end

function Base.full{T}(B::BrokenArrowBidiagonal{T})
    n = size(B, 1)
    k = length(B.av)
    M = zeros(T, n, n)
    for i=1:n
        M[i, i] = B.dv[i]
    end
    for i=1:k
        M[i,k+1] = B.av[i]
    end
    for i=k+1:n-1
        M[i, i+1] = B.ev[i-k]
    end
    M
end

Base.svdfact(B::BrokenArrowBidiagonal) = svdfact(full(B)) #XXX This can be much faster

"""
    PartialFactorization

Partial factorization object which is an approximation of a matrix.

    A ≈ P * [B 0; 0 β] * Q

"""
type PartialFactorization{T, Tr} <: Factorization{T}
    P :: Matrix{T}
    Q :: Matrix{T}
    B :: AbstractMatrix{Tr}
    β :: Tr
end

####################
# API method calls #
####################

"""
    svdl(Master, A)

Compute some singular values (and optionally vectors) using Golub-Kahan-Lanczos
bidiagonalization \cite{Golub1965} with thick restarting \cite{Wu2000}.

# Arguents

* `A` : The matrix or matrix-like object whose singular values are desired.

## Keywords

* `nsv::Int = 6`: number of singular values requested.

* `v0 = random unit vector`: starting guess vector in the domain of `A`.
The length of `q` should be the number of columns in `A`.

* `k::Int = 2nsv`: maximum number of Lanczos vectors to compute before restarting.

* `j::Int = nsv`: number of vectors to keep at the end of the restart.
We don't recommend j < nsv.

* `maxiter::Int = minimum(size(A))`: maximum number of iterations to run.

* `verbose::Bool = false`: print information at each iteration.

* `tol::Real = √eps()`: maximum absolute error in each desired singular value.

* `reltol::Real=√eps()`: maximum error in each desired singular value relative to the
estimated norm of the input matrix.

* `method::Symbol=:ritz`: restarting algorithm to use. Valid choices are:
    - `:ritz`: Thick restart with Ritz values [Wu2000].
    - `:harmonic`: Restart with harmonic Ritz values [Baglama2005].

* `vecs::Symbol = :none`: singular vectors to return.
    - `:both`: Both left and right singular vectors are returned.
    - `:left`: Only the left singular vectors are returned.
    - `:right`: Only the right singular vectors are returned.
    - `:none`: No singular vectors are returned.

* `dolock::Bool=false`: If `true`, locks converged Ritz values, removing them
from the Krylov subspace being searched in the next macroiteration.

# Output

* `Σ`: list of the desired singular values if `vecs == :none` (the default),
    otherwise returns an `SVD` object with the desired singular vectors filled in.

* `L`: computed partial factorizations of A.

* `ch::ConvergenceHistory`:
    - `:betas` => `betas`: The history of the computed betas.
    - `:Bs` => `Bs`: The history of the computed projected matrices.
    - `:ritz` => `ritzvalhist`: Ritz values computed at each iteration.
    - `:conv` => `convhist`: Convergence data.

# Implementation notes

The implementation of thick restarting follows closely that of SLEPc as
described in [Hernandez2008]. Thick restarting can be turned off by setting `k
= maxiter`, but most of the time this is not desirable.

The singular vectors are computed directly by forming the Ritz vectors from the
product of the Lanczos vectors `L.P`/`L.Q` and the singular vectors of `L.B`.
Additional accuracy in the singular triples can be obtained using inverse
iteration.

# References

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

@article{Wu2000,
    author = {Wu, Kesheng and Simon, Horst},
    journal = {SIAM Journal on Matrix Analysis and Applications},
    number = 2,
    pages = {602--616},
    title = {Thick-Restart {L}anczos Method for Large Symmetric Eigenvalue Problems},
    volume = 22,
    year = 2000
}

@article{Baglama2005,
    author = {Baglama, James and Reichel, Lothar},
    doi = {10.1137/04060593X},
    journal = {SIAM Journal on Scientific Computing},
    number = 1,
    pages = {19--42},
    title = {Augmented Implicitly Restarted {L}anczos Bidiagonalization Methods},
    volume = 27,
    year = 2005
}

@article{Hernandez2008,
    author = {Hern\'{a}ndez, Vicente and Rom\'{a}n, Jos\'{e} E and Tom\'{a}s,
    Andr\'{e}s},
    journal = {Electronic Transactions on Numerical Analysis},
    pages = {68--85},
    title = {A Robust and Efficient Parallel {SVD} Solver based on Restarted
        {L}anczos Bidiagonalization},
    url = {http://etna.mcs.kent.edu/volumes/2001-2010/vol31/abstract.php?vol=31\&pages=68-85},
    volume = 31,
    year = 2008
}
```

"""
svdl(A; kwargs...) = svdl_method(A; kwargs...)

function svdl(::Type{Master}, A;
    tol::Real=√eps(), plot::Bool=false, nsv::Int=6, k::Int=2nsv,
    maxiter::Int=minimum(size(A)), method::Symbol=:ritz, kwargs...
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:conv, maxiter, T=BitArray)
    reserve!(log,:ritz, maxiter, k)
    reserve!(log,:resnorm, maxiter, k)
    Bs_type = (method == :ritz) ? BrokenArrowBidiagonal : UpperTriangular
    reserve!(log,:Bs, maxiter, T=Bs_type)
    reserve!(log,:betas, maxiter)
    X, L = svdl_method(A;
        tol=tol, log=log, k=k, nsv=nsv, maxiter=maxiter, method=method, kwargs...)
    shrink!(log)
    plot && showplot(log)
    X, L, log
end

#########################
# Method Implementation #
#########################

function svdl_method(A;
    nsv::Int=6, k::Int=2nsv, j::Int=nsv,
    v0::AbstractVector = Vector{eltype(A)}(randn(size(A, 2))) |> x->scale!(x, inv(norm(x))),
    maxiter::Int=minimum(size(A)), tol::Real=√eps(), reltol::Real=√eps(),
    verbose::Bool=false, method::Symbol=:ritz, vecs=:none,
    dolock::Bool=false, log::MethodLog=DummyHistory()
    )

    T0 = time_ns()
    @assert k>1
    L = build(A, v0, k)

    converged = false

    iter = 0
    local F
    for iter in 1:maxiter
        #@assert size(L.B) == (k, k)
        F = svdfact(L.B)::LinAlg.SVD{eltype(v0), typeof(real(one(eltype(v0)))), Matrix{eltype(v0)}}
        iter==1 && @assert eltype(F)==eltype(v0)
        if method == :ritz
            thickrestart!(A, L, F, j)
        elseif method == :harmonic
            harmonicrestart!(A, L, F, j)
        else
            throw(ArgumentError("Unknown restart method $method"))
        end
        extend!(A, L, k)
        if verbose
	        elapsedtime = round((time_ns()-T0)*1e-9, 3)
	        info("Iteration $iter: $elapsedtime seconds")
        end

        nextiter!(log)
        conv = isconverged(L, F, nsv, tol, reltol, log, verbose)

        push!(log, :conv, conv)
        push!(log, :ritz, F[:S][1:k])
        push!(log, :Bs, deepcopy(L.B))
        push!(log, :betas, L.β)

        #Lock
        if method == :ritz && dolock
            for i in eachindex(conv)
                if conv[i]
                    L.B.av[i] = 0
                end
            end
        end
        all(conv) && (setconv(log, true); break)
    end
    setmvps(log, iter)

    #Compute singular vectors as necessary and return them in the output
    values = F[:S][1:nsv]
    m, n = size(A)

    leftvecs = if vecs == :left || vecs == :both
        L.P*sub(F[:U], :, 1:nsv)
    else
        zeros(eltype(v0), m, 0)
    end

    rightvecs = if vecs == :right || vecs == :both
        (sub(L.Q, :, 1:size(L.Q,2)-1)*sub(F[:V], :, 1:nsv))'
    else
        zeros(eltype(v0), 0, n)
    end

    if vecs == :none
        values, L, converged
    else
        LinAlg.SVD(leftvecs, values, rightvecs), L, converged
    end
end

"""
    isconverged(L, F, k, tol, reltol, log, verbose=false)

Determine if any singular values in a partial factorization have converged.

# Arguments

* `L::PartialFactorization` : a `PartialFactorization` computed by an iterative
method such as `svdl`.

* `F::Base.LinAlg.SVD`: a `SVD` factorization computed for `L.B`.

* `k::Int` : number of singular values to check.

* `tol::Real`: absolute tolerance for a Ritz value to be considered converged.

* `reltol::Real`: relative tolerance for a Ritz value to be considered converged.

* `verbose::Bool = false`: if `true`, prints out all the results of convergence tests.

# Implementation notes

This convergence test routine uses a variety of different tests.

1. The crudest estimate of the error bound is a simple error bound which dates
back at least to Wilkinson's classic book [Wilkinson1965:Ch.3 §53 p.170]

```bibtex
@book{Wilkinson1965,
    address = {Oxford, UK},
    author = {J H Wilkinson},
    publisher = {Oxford},
    title = {The Algebraic Eigenvalue Problem},
    year = 1965
}
```

2. When the Ritz values become sufficiently well-separated, more refined estimates
can be derived from the Rayleigh-Ritz properties of the Krylov process, as
described in [Wilkinson1965:Ch.3 §54-55 p.173, Yamamoto1980, Ortega1990]

```bibtex
@book{Ortega1990,
    address = {Philadelphia, PA},
    author = {Ortega, James M},
    doi = {10.1137/1.9781611971323},
    edition = {2},
    publisher = {SIAM},
    series = {Classics in Applied Mathematics},
    title = {Numerical Analysis: A Second Course},
    url = {http://epubs.siam.org/doi/book/10.1137/1.9781611971323},
    year = 1990
}

@article{Yamamoto1980,
    author = {Yamamoto, Tetsuro},
    doi = {10.1007/BF01396059},
    journal = {Numerische Mathematik},
    number = 2,
    pages = {189--199},
    title = {Error bounds for computed eigenvalues and eigenvectors},
    volume = {34},
    year = 1980
}

@article{Geurts1982,
    author = {Geurts, A J},
    doi = {10.1007/BF01399313},
    journal = {Numerische Mathematik},
    month = feb,
    number = {1},
    pages = {85--96},
    title = {A contribution to the theory of condition},
    volume = {39},
    year = {1982}
}

@article{Deif1989,
    author = {Deif, A.},
    doi = {10.1007/BF01396348},
    journal = {Numerische Mathematik},
    month = jun,
    number = {6},
    pages = {625--626},
    title = {A relative backward perturbation theorem for the eigenvalue problem},
    volume = {56},
    year = {1989}
}
```
"""
function isconverged(L::PartialFactorization,
        F::Base.LinAlg.SVD, k::Int, tol::Real, reltol::Real,
        log::MethodLog, verbose::Bool=false)

    @assert tol ≥ 0

    σ = F[:S][1:k]
    Δσ= L.β*abs(F[:U][end, 1:k])

    #Best available eigenvalue bounds
    δσ = copy(Δσ)

    #Update better error bounds from Rayleigh-Ritz
    #In theory these only apply if we know all the values
    #But so long as the gap between the converged Ritz values and all the
    #others is larger than the smallest empirical spectral gap d then we're fine.
    #Can also get some eigenvector statistics!!
    if k > 1
        d = Inf
        for i in eachindex(σ), j=1:i-1
            d = min(d, abs(σ[i] - σ[j]))
        end
        verbose && println("Smallest empirical spectral gap: ", d)
        #Chatelein 1993 - normwise backward error associated with approximate invariant subspace
        #Use the largest singular (Ritz) value to estimate the 2-norm of the matrix
        verbose && println("Normwise backward error associated with subspace: ", L.β/σ[1])
        for i in eachindex(Δσ)
            α = Δσ[i]

            #Simple error bound
            if 2α ≤ d
                #Rayleigh-Ritz bounds
                x = abs(α/(d-α)*√(1+(α/(d-α))^2))
                verbose && println("Rayleigh-Ritz error bound on eigenvector: $x")
                #δσ[i] = min(δσ[i], x)

                y = α^2/d #[Wilkinson:Ch.3 Appendix (4), p.188]
                verbose && println("Rayleigh-Ritz error bound on eigenvalue: $y")
                δσ[i] = min(δσ[i], y)
            end

	    verbose && println("Ritz value ", i, ": ", σ[i], " ± ", signif(δσ[i], 3))

            #Estimate of the normwise backward error [Deif 1989]
            ##Use the largest singular (Ritz) value to estimate the 2-norm of the matrix
            verbose && println("Normwise backward error estimate: ", α/σ[1])

            #Geurts, 1982 - componentwise backward error also known
            #A. J. Geurts, (1982), A contribution to the theory of condition,
            #Numer.  Math., #39, 85-96.
        end
    end

    #Estimate condition number and see if two-sided reorthog is needed
    if verbose && (F[:S][1]/F[:S][end]) > 1/√eps()
        warn("Two-sided reorthogonalization should be used but is not implemented")
    end

    push!(log, :resnorm, δσ[1:k])
    conv = (δσ[1:k] .< max(tol, reltol*σ[1]))::BitVector
end

"""
    build(A, q, k)

# References

\cite{Hernandez2008}

"""
function build{T}(A, q::AbstractVector{T}, k::Int)
    m, n = size(A)
    Tr = typeof(real(one(T)))
    β = norm(q)
    scale!(q, inv(β))
    p = A*q
    α = convert(Tr, norm(p))
    scale!(p, inv(α))
    extend!(A, PartialFactorization(
        reshape(p, m, 1), reshape(q, n, 1), Bidiagonal([α], Tr[], true), β
        ), k)
end


"""
    thickrestart!(A, L, F, l)

Thick restart (with ordinary Ritz values).

# Reference

[Hernandez2008]

"""
function thickrestart!{T,Tr}(A, L::PartialFactorization{T,Tr},
        F::Base.LinAlg.SVD{Tr,Tr}, l::Int)

    k = size(F[:V], 1)
    m, n = size(A)
    #@assert size(L.P) == (m, k)
    #@assert size(L.Q) == (n, k+1)

    Q = sub(L.Q, :,1:k)*sub(F[:V], :,1:l)
    L.Q = [Q sub(L.Q, :, k+1)]
    #Be pedantic about ensuring normalization
    #L.Q = qr(L.Q)[1]
    #@assert all([norm(L.Q[:,i]) ≈ 1 for i=1:size(L.Q,2)])

    f = A*sub(L.Q, :, l+1)
    ρ = L.β * reshape(F[:U][end, 1:l], l)
    L.P = sub(L.P, :, 1:k)*sub(F[:U], :, 1:l)

    #@assert ρ[i] ≈ f⋅L.P[:, i]
    f -= L.P*ρ
    α = convert(Tr, norm(f))
    scale!(f, inv(α))
    L.P = [L.P f]

    g = A'f - α*L.Q[:, end]
    L.β = β = norm(g)
    L.B = BrokenArrowBidiagonal([F[:S][1:l]; α], ρ, typeof(β)[])
    #@assert size(L.P, 2) == size(L.Q, 2) == size(L.B, 2)
    L
end

"""
    harmonicrestart!(A, L, F, k)

Thick restart with harmonic Ritz values.

# Reference

[Baglama2005] - note that they have P and Q swapped relative to our notation,
which follows that of [Hernandez2008]

"""
function harmonicrestart!{T,Tr}(A, L::PartialFactorization{T,Tr},
        F::Base.LinAlg.SVD{Tr,Tr}, k::Int)

    m = size(L.B, 1)::Int
    #@assert size(L.P,2)==m==size(L.Q,2)-1

    F0 = F# svdfact(L.B)
    ρ = L.β*F0[:U][end,:] #Residuals of singular values


    #Construct broken arrow matrix
    if VERSION < v"0.5.0-"
        BA = [diagm(F0[:S]) ρ']
    else #slicing behavior changed in 0.5
        BA = [diagm(F0[:S]) ρ]
    end
    F2 = svdfact!(BA, thin=false)

    #Take k largest triplets
    Σ = (F2[:S]::Vector{Tr})[1:k]
    U = F0[:U]*sub(F2[:U],:,1:k)
    M = eye(T, m+1)
    M[1:m, 1:m] = F0[:V]::Matrix{T}
    M = M * F2[:V]
    Mend = M[end, 1:k]
    #Compute scaled residual from the harmonic Ritz problem
    r0 = zeros(Tr, m)
    r0[end] = 1
    isa(L.B, Bidiagonal) && @assert L.B.isupper
    r = try
        #(L.B\r0)
        A_ldiv_B!(L.B, r0)
    catch exc
        if isa(exc, LinAlg.LAPACKException) ||
                isa(exc, LinAlg.SingularException) #B\r is singular

            pinv(full(L.B))*r0
        else rethrow(exc) end
    end::Vector{Tr}
    scale!(r, L.β)
    M::Matrix{T} = sub(M,1:m, :) + r*sub(M,m+1,:)

    M2 = zeros(T, m+1, k+1)
    M2[1:m, 1:k] = M[:,1:k]
    M2[1:m, k+1] = -r
    M2[m+1, k+1] = 1
    Q, R = qr(M2)

    Q = L.Q*Q
    P = L.P*sub(U,:,1:k)

    if VERSION < v"0.5.0-"
        R = sub(R,1:k+1,1:k) + sub(R,:,k+1)*Mend
    else
        R = sub(R,1:k+1,1:k) + sub(R,:,k+1)*Mend'
    end

    f = A*sub(Q,:,k+1)
    f -= P*(P'f)
    α = convert(Tr, norm(f))
    scale!(f, inv(α))
    P = [P f]
    B = UpperTriangular{Tr,Matrix{Tr}}([(Diagonal(Σ)*triu(R')); zeros(Tr,1,k) α])
    #@assert size(P, 2) == size(B, 1) == size(Q, 2)
    g = A'f
    q = sub(Q,:,k+1)
    #g-= (g⋅q)*q
    axpy!(-(g⋅q), q, g)
    β = convert(Tr, norm(g))
    scale!(g, inv(β))
    #@assert size(P, 2) == size(Q, 2) == size(B, 2)
    L.P = P
    L.Q = Q
    L.B = B
    L.β = β
    L
end

"""
    extend!{T,Tr}(A, L, k orthleft, orthright, α)

Extend a PartialFactorization L using GKL bidiagonalization with k extra pairs
of Lanczos vectors.

# Arguments

* `A`: matrix or linear map generating the Lanczos vectors

* `L::PartialFactorization`: partial factorization.

* `orthleft::Bool = false`: orthogonalize left Lanczos vectors.

* `orthright::Bool = true`: orthogonalize right Lanczos vectors.

* `α::Real = 1/√2`: criterion for doing a second reorthogonalization.

# Implementation notes

The implementation mostly follows the description in [Simon2000,Hernandez2008].

The reorthogonalization method used is using double classical Gram-Schmidt
full reorthogonalization. As explained in the numerical analysis literature by
Kahan, Golub, Rutishauser, and others in the 1970s, double classical
Gram-Schmidt reorthogonalization always suffices to keep vectors orthogonal to
within machine precision. As described in [Bjorck2015], α is a threshold
for determinining when the second orthogonalization is necessary. -log10(α) is
the number of (decimal) digits lost due to cancellation. Common choices are
α=0.1 [Rutishauser] and α=1/√2 [Daniel1976] (our default).

In most situations it suffices to orthogonalize either the left vectors or the
right vectors, except when the matrix norm exceeds `1/√eps(eltype(A))`, in
which case it will be necessary to orthogonalize both sets of vectors. See
[Simon2000].

```bibtex
@book{Bjorck2015,
    author = {Bj{\"{o}}rck, {\AA}ke},
    doi = {10.1007/978-3-319-05089-8},
    publisher = {Springer},
    series = {Texts in Applied Mathematics},
    title = {Numerical Methods in Matrix Computations},
    year = {2015}
}

@article{Simon2000,
    author = {Simon, Horst D. and Zha, Hongyuan},
    doi = {10.1137/S1064827597327309},
    journal = {SIAM Journal on Scientific Computing},
    number = {6},
    pages = {2257--2274},
    title = {Low-Rank Matrix Approximation Using the {Lanczos} Bidiagonalization Process with Applications},
    volume = {21},
    year = {2000}
}

@article{Daniel1976,
    author = {Daniel, J. W. and Gragg, W. B. and Kaufman, L. and Stewart, G. W.},
    doi = {10.1090/S0025-5718-1976-0431641-8},
    journal = {Mathematics of Computation},
    number = {136},
    pages = {772--795},
    title = {Reorthogonalization and stable algorithms for updating the {Gram-Schmidt QR} factorization},
    volume = {30},
    year = {1976}
}
```

"""
function extend!{T,Tr}(A, L::PartialFactorization{T, Tr}, k::Int,
    orthleft::Bool=false, orthright::Bool=true, α::Real = 1/√2)

    l = size(L.B, 2)::Int-1
    p = L.P[:,l+1]

    m, n = size(A)
    q = zeros(T, n)
    #@assert l+1 == size(L.B, 1) == size(L.B, 2) == size(L.P, 2)

    if isa(L.B, UpperTriangular) #Cannot be appended
        B = zeros(Tr, k, k)
        B[1:size(L.B,1), 1:size(L.B,2)] = L.B
        L.B = UpperTriangular(B)
        #@assert size(L.B) == (k, k)
    end

    β = L.β
    for j=l+1:k
        Ac_mul_B!(q, A, p) #q = A'p

        if orthright #Orthogonalize right Lanczos vector
            #Do double classical Gram-Schmidt reorthogonalization
            oldqnorm = norm(q)
            q -= L.Q*(L.Q'q)
            if norm(q) ≤ α * oldqnorm
                q -= L.Q*(L.Q'q)
            end
        end

        β = norm(q)
        scale!(q, inv(β))

        L.Q = [L.Q q]
        j==k && break

        #p = A*q - β*p
        A_mul_B!(p, A, q)
        Base.LinAlg.axpy!(-β, sub(L.P, :, j), p)

        if orthleft #Orthogonalize left Lanczos vector
            #Do double classical Gram-Schmidt reorthogonalization
            oldpnorm = norm(p)
            p -= L.P*(L.P'p)
            if norm(p) ≤ α * oldpnorm
                p -= L.P*(L.P'p)
            end
        end

        α = norm(p)
        scale!(p, inv(α))
        if isa(L.B, Bidiagonal) || isa(L.B, BrokenArrowBidiagonal)
            push!(L.B.dv, α)
            push!(L.B.ev, β)
        else
            L.B[j+1, j+1] = α
            L.B[j  , j+1] = β
        end
        L.P = [L.P p]
    end
    L.β = β
    L
end

let
    B = BrokenArrowBidiagonal([1, 2, 3], [1, 2], Int[])
    @assert full(B) == [1 0 1; 0 2 2; 0 0 3]
end
