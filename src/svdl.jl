export svdl

import Base: size, getindex, Matrix
import LinearAlgebra.svd
using LinearAlgebra

"""
Matrix of the form
```
 d1          a_1
    d2       a_2
      ...    ...
         d_l a_l
             d_l+1 e_l+1
                   ...  e_k-1
                        d_k
```
"""
mutable struct BrokenArrowBidiagonal{T} <: AbstractMatrix{T}
    dv::Vector{T}
    av::Vector{T}
    ev::Vector{T}
end

size(B::BrokenArrowBidiagonal) = (n=length(B.dv); (n, n))

function size(B::BrokenArrowBidiagonal, n::Int)
    if n==1 || n==2
        return length(B.dv)
    else
        throw(ArgumentError("invalid dimension $n"))
    end
end

function getindex(B::BrokenArrowBidiagonal{T}, i::Int, j::Int) where {T}
    n = size(B, 1)
    k = length(B.av)
    if !(1 ≤ i ≤ n && 1 ≤ j ≤ n)
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

function Matrix(B::BrokenArrowBidiagonal{T}) where {T}
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

svd(B::BrokenArrowBidiagonal) = svd(Matrix(B)) #XXX This can be much faster

"""
Partial factorization object which is an approximation of a matrix

    A ≈ P * [B 0; 0 β] * Q
"""
mutable struct PartialFactorization{T, Tr} <: Factorization{T}
    P :: Matrix{T}
    Q :: Matrix{T}
    B :: AbstractMatrix{Tr}
    β :: Tr
end

####################
# API method calls #
####################

"""
    svdl(A) -> Σ, L, [history]

Compute some singular values (and optionally vectors) using Golub-Kahan-Lanczos
bidiagonalization [^Golub1965] with thick restarting [^Wu2000].

If `log` is set to `true` is given, method will output a tuple `X, L, ch`. Where
`ch` is a `ConvergenceHistory` object. Otherwise it will only return `X, L`.

# Arguments

- `A` : The matrix or matrix-like object whose singular values are desired.

## Keywords

- `nsv::Int = 6`: number of singular values requested;
- `v0 = random unit vector`: starting guess vector in the domain of `A`.
  The length of `q` should be the number of columns in `A`;
- `k::Int = 2nsv`: maximum number of Lanczos vectors to compute before restarting;
- `j::Int = nsv`: number of vectors to keep at the end of the restart.
  We don't recommend j < nsv;
- `maxiter::Int = minimum(size(A))`: maximum number of iterations to run;
- `verbose::Bool = false`: print information at each iteration;
- `tol::Real = √eps()`: maximum absolute error in each desired singular value;
- `reltol::Real=√eps()`: maximum error in each desired singular value relative to the
  estimated norm of the input matrix;
- `method::Symbol=:ritz`: restarting algorithm to use. Valid choices are:
  1. `:ritz`: Thick restart with Ritz values [Wu2000].
  2. `:harmonic`: Restart with harmonic Ritz values [Baglama2005].
- `vecs::Symbol = :none`: singular vectors to return.
  1. `:both`: Both left and right singular vectors are returned.
  2. `:left`: Only the left singular vectors are returned.
  3. `:right`: Only the right singular vectors are returned.
  4. `:none`: No singular vectors are returned.
- `dolock::Bool=false`: If `true`, locks converged Ritz values, removing them
  from the Krylov subspace being searched in the next macroiteration;
- `log::Bool = false`: output an extra element of type `ConvergenceHistory`
  containing extra information of the method execution.

# Return values

**if `log` is `false`**

- `Σ`: list of the desired singular values if `vecs == :none` (the default), otherwise
  returns an `SVD` object with the desired singular vectors filled in;
- `L`: computed partial factorizations of A.

**if `log` is `true`**

- `Σ`: list of the desired singular values if `vecs == :none` (the default),
otherwise returns an `SVD` object with the desired singular vectors filled in;
- `L`: computed partial factorizations of A;
- `history`: convergence history.

**ConvergenceHistory keys**

- `:betas` => `betas`: The history of the computed betas.
- `:Bs` => `Bs`: The history of the computed projected matrices.
- `:ritz` => `ritzvalhist`: Ritz values computed at each iteration.
- `:conv` => `convhist`: Convergence data.

[^Golub1965]:
    Golub, Gene, and William Kahan. "Calculating the singular values and pseudo-inverse
    of a matrix." Journal of the Society for Industrial and Applied Mathematics,
    Series B: Numerical Analysis 2.2 (1965): 205-224.
[^Wu2000]:
    Wu, Kesheng, and Horst Simon. "Thick-restart Lanczos method for large symmetric
    eigenvalue problems." SIAM Journal on Matrix Analysis and Applications 22.2
    (2000): 602-616.
"""
function svdl(A;
    nsv::Int=6, k::Int=2nsv, tol::Real=√eps(),
    maxiter::Int=minimum(size(A)), method::Symbol=:ritz, log::Bool=false, kwargs...
    )
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(BitArray, history,:conv, maxiter)
    reserve!(history,[:ritz,:resnorm], maxiter, k)
    Bs_type = (method == :ritz) ? BrokenArrowBidiagonal : UpperTriangular
    reserve!(Bs_type, history,:Bs, maxiter)
    reserve!(history,:betas, maxiter)
    X, L = svdl_method!(history, A, nsv; k=k, tol=tol, maxiter=maxiter, method=method, kwargs...)
    log && shrink!(history)
    log ? (X, L, history) : (X, L)
end

#########################
# Method Implementation #
#########################

function svdl_method!(log::ConvergenceHistory, A, l::Int=min(6, size(A,1)); k::Int=2l,
    j::Int=l, v0::AbstractVector = Vector{eltype(A)}(randn(size(A, 2))) |> x->rmul!(x, inv(norm(x))),
    maxiter::Int=minimum(size(A)), tol::Real=√eps(), reltol::Real=√eps(),
    verbose::Bool=false, method::Symbol=:ritz, vecs=:none, dolock::Bool=false)

    T0 = time_ns()
    @assert k>1
    L = build(log, A, v0, k)

    iter = 0
    local F
    for iter in 1:maxiter
        nextiter!(log)

        #@assert size(L.B) == (k, k)
        F = svd(L.B)::LinearAlgebra.SVD{eltype(v0), typeof(real(one(eltype(v0)))), Matrix{eltype(v0)}}
        iter==1 && @assert eltype(F)==eltype(v0)
        if method == :ritz
            thickrestart!(A, L, F, j)
        elseif method == :harmonic
            harmonicrestart!(A, L, F, j)
        else
            throw(ArgumentError("Unknown restart method $method"))
        end
        extend!(log, A, L, k)
        if verbose
            elapsedtime = round((time_ns()-T0)*1e-9, digits=3)
            info("Iteration $iter: $elapsedtime seconds")
        end

        conv = isconverged(L, F, l, tol, reltol, log, verbose)

        push!(log, :conv, conv)
        push!(log, :ritz, F.S[1:k])
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
    shrink!(log)

    #Compute singular vectors as necessary and return them in the output
    values = F.S[1:l]
    m, n = size(A)

    leftvecs = if vecs == :left || vecs == :both
        L.P*view(F.U, :, 1:l)
    else
        zeros(eltype(v0), m, 0)
    end

    rightvecs = if vecs == :right || vecs == :both
        (view(L.Q, :, 1:size(L.Q,2)-1)*view(F.V, :, 1:l))'
    else
        zeros(eltype(v0), 0, n)
    end

    if vecs == :none
        values, L
    else
        LinearAlgebra.SVD(leftvecs, values, rightvecs), L
    end
end

"""
    isconverged(L, F, k, tol, reltol, log, verbose=false)

Determine if any singular values in a partial factorization have converged.

# Arguments

- `L::PartialFactorization`: a `PartialFactorization` computed by an iterative
method such as `svdl`;
- `F::LinearAlgebra.SVD`: a `SVD` factorization computed for `L.B`;
- `k::Int` : number of singular values to check;
- `tol::Real`: absolute tolerance for a Ritz value to be considered converged;
- `reltol::Real`: relative tolerance for a Ritz value to be considered converged;
- `verbose::Bool = false`: if `true`, prints out all the results of convergence tests.

# Implementation note

This convergence test routine uses a variety of different tests.
1. The crudest estimate of the error bound is a simple error bound which dates
back at least to Wilkinson's classic book [^Wilkinson1965] Ch.3 §53 p.170.
2. When the Ritz values become sufficiently well-separated, more refined estimates
can be derived from the Rayleigh-Ritz properties of the Krylov process, as
described in [^Wilkinson1965] Ch.3 §54-55 p.173, [^Yamamoto1980], [^Ortega1990],
[^Geurts1982] and [^Deif1989].

[^Wilkinson1965]:
    Wilkinson, James Hardy. The algebraic eigenvalue problem.
    Vol. 87. Oxford: Clarendon Press, 1965.
[^Yamamoto1980]:
    Yamamoto, Tetsuro. "Error bounds for computed eigenvalues
    and eigenvectors." Numerische Mathematik 34.2 (1980): 189-199.
[^Ortega1990]:
    Ortega, James M. Numerical analysis: a second course.
    Society for Industrial and Applied Mathematics, 1990.
[^Geurts1982]:
    Geurts, A. J. "A contribution to the theory of condition."
    Numerische Mathematik 39.1 (1982): 85-96.
[^Deif1989]:
    Deif, A. "A relative backward perturbation theorem for the eigenvalue
    problem." Numerische Mathematik 56.6 (1989): 625-626.
"""
function isconverged(L::PartialFactorization, F::LinearAlgebra.SVD, k::Int, tol::Real,
    reltol::Real, log::ConvergenceHistory, verbose::Bool=false
    )

    @assert tol ≥ 0

    σ = F.S[1:k]
    Δσ= L.β * abs.(F.U[end, 1 : k])

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

	    verbose && println("Ritz value ", i, ": ", σ[i], " ± ", round(δσ[i], sigdigits=3))

            #Estimate of the normwise backward error [Deif 1989]
            ##Use the largest singular (Ritz) value to estimate the 2-norm of the matrix
            verbose && println("Normwise backward error estimate: ", α/σ[1])

            #Geurts, 1982 - componentwise backward error also known
            #A. J. Geurts, (1982), A contribution to the theory of condition,
            #Numer.  Math., #39, 85-96.
        end
    end

    #Estimate condition number and see if two-sided reorthog is needed
    if verbose && (F.S[1]/F.S[end]) > 1/√eps()
        @warn("Two-sided reorthogonalization should be used but is not implemented")
    end

    push!(log, :resnorm, δσ[1:k])
    conv = (δσ[1:k] .< max(tol, reltol*σ[1]))::BitVector
end

#Hernandez2008
function build(log::ConvergenceHistory, A, q::AbstractVector{T}, k::Int) where {T}
    m, n = size(A)
    Tr = typeof(real(one(T)))
    β = norm(q)
    q .*= inv(β)
    p = A*q
    α = convert(Tr, norm(p))
    p .*= inv(α)
    bidiag = Bidiagonal([α], Tr[], :U)
    extend!(log, A, PartialFactorization(reshape(p, m, 1), reshape(q, n, 1), bidiag, β), k)
end


"""
    thickrestart!(A, L,F, l)

Thick restart (with ordinary Ritz values)

# References

[^Hernandez2008]

"""
function thickrestart!(A, L::PartialFactorization{T,Tr},
        F::LinearAlgebra.SVD{Tr,Tr}, l::Int) where {T,Tr}

    k = size(F.V, 1)
    m, n = size(A)
    #@assert size(L.P) == (m, k)
    #@assert size(L.Q) == (n, k+1)

    Q = view(L.Q, :,1:k)*view(F.V, :,1:l)
    L.Q = [Q view(L.Q, :, k+1)]
    #Be pedantic about ensuring normalization
    #L.Q = qr(L.Q)[1]
    #@assert all([norm(L.Q[:,i]) ≈ 1 for i=1:size(L.Q,2)])

    f = A*view(L.Q, :, l+1)
    ρ = L.β * reshape(F.U[end, 1:l], l)
    L.P = view(L.P, :, 1:k)*view(F.U, :, 1:l)

    #@assert ρ[i] ≈ f⋅L.P[:, i]
    f -= L.P*ρ
    α = convert(Tr, norm(f))
    f .*= inv(α)
    L.P = [L.P f]

    g = A'f - α*L.Q[:, end]
    L.β = β = norm(g)
    L.B = BrokenArrowBidiagonal([F.S[1:l]; α], ρ, typeof(β)[])
    #@assert size(L.P, 2) == size(L.Q, 2) == size(L.B, 2)
    L
end

"""
    harmonicrestart!(A, L, F, k) -> L

Thick restart with harmonic Ritz values. See [^Baglama2005] but note that
they have P and Q swapped relative to our notation, which follows that
of [^Hernandez2008]

[^Baglama2005]:
    Baglama, James, and Lothar Reichel. "Augmented implicitly restarted
    Lanczos bidiagonalization methods." SIAM Journal on Scientific
    Computing 27.1 (2005): 19-42.
[^Hernandez2008]:
    Vicente Hernández, José E. Román, and Andrés Tomás. "A robust and
    efficient parallel SVD solver based on restarted Lanczos bidiagonalization."
    Electronic Transactions on Numerical Analysis 31 (2008): 68-85.

"""
function harmonicrestart!(A, L::PartialFactorization{T,Tr},
        F::LinearAlgebra.SVD{Tr,Tr}, k::Int) where {T,Tr}

    m = size(L.B, 1)::Int
    #@assert size(L.P,2)==m==size(L.Q,2)-1

    F0 = F# svd(L.B)
    ρ = L.β*F0.U[end,:] #Residuals of singular values


    #Construct broken arrow matrix
    BA = [Matrix(Diagonal(F0.S)) ρ]
    F2 = svd!(BA; full=true)

    #Take k largest triplets
    Σ = (F2.S::Vector{Tr})[1:k]
    U = F0.U*view(F2.U,:,1:k)
    M = Matrix{T}(I, m+1, m+1)
    M[1:m, 1:m] = F0.V::Union{Matrix{T}, Adjoint{T,Matrix{T}}}
    M = M * F2.V
    Mend = M[end, 1:k]
    #Compute scaled residual from the harmonic Ritz problem
    r0 = zeros(Tr, m)
    r0[end] = 1
    if isa(L.B, Bidiagonal)
        @assert L.B.uplo == 'U'
    end
    r = try
        #(L.B\r0)
        ldiv!(L.B, r0)
    catch exc
        if isa(exc, LinearAlgebra.LAPACKException) ||
            isa(exc, LinearAlgebra.SingularException) #B\r is singular
            pinv(Matrix(L.B))*r0
        else rethrow(exc) end
    end::Vector{Tr}
    r .*= L.β
    M::Matrix{T} = view(M,1:m, :) + r*view(M,m+1:m+1,:)

    M2 = zeros(T, m+1, k+1)
    M2[1:m, 1:k] = M[:,1:k]
    M2[1:m, k+1] = -r
    M2[m+1, k+1] = 1
    QRF = qr(M2)
    Q, R = QRF.Q, QRF.R

    Q = L.Q*view(Q,:,1:k+1)
    P = L.P*view(U,:,1:k)

    R = view(R,1:k+1,1:k) + view(R,:,k+1)*Mend'

    f = A*view(Q,:,k+1)
    f -= P*(P'f)
    α = convert(Tr, norm(f))
    f .*= inv(α)
    P = [P f]
    B = UpperTriangular{Tr,Matrix{Tr}}([(Diagonal(Σ)*triu(R')); zeros(Tr,1,k) α])
    #@assert size(P, 2) == size(B, 1) == size(Q, 2)
    g = A'f
    q = view(Q,:,k+1)
    #g-= (g⋅q)*q
    g .-= (g⋅q)*q
    β = convert(Tr, norm(g))
    g .*= inv(β)
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

- `A`: matrix or linear map generating the Lanczos vectors;
- `L::PartialFactorization`: partial factorization;
- `orthleft::Bool = false`: orthogonalize left Lanczos vectors;
- `orthright::Bool = true`: orthogonalize right Lanczos vectors;
- `α::Real = 1/√2`: criterion for doing a second reorthogonalization.

## Implementation details

The implementation mostly follows the description in [^Simon2000] and [^Hernandez2008].

The reorthogonalization method used is using double classical Gram-Schmidt
full reorthogonalization. As explained in the numerical analysis literature by
Kahan, Golub, Rutishauser, and others in the 1970s, double classical
Gram-Schmidt reorthogonalization always suffices to keep vectors orthogonal to
within machine precision. As described in [^Bjorck2015], `α` is a threshold
for determinining when the second orthogonalization is necessary. -log10(α) is
the number of (decimal) digits lost due to cancellation. Common choices are
`α=0.1` [Rutishauser] and `α=1/√2` [Daniel1976] (our default).

In most situations it suffices to orthogonalize either the left vectors or the
right vectors, except when the matrix norm exceeds `1/√eps(eltype(A))`, in
which case it will be necessary to orthogonalize both sets of vectors. See
[^Simon2000].

[^Bjorck2015]:
    Björck, Åke. Numerical methods in matrix computations.
    New York: Springer, 2015.
[^Simon2000]:
    Simon, Horst D., and Hongyuan Zha. "Low-rank matrix approximation
    using the Lanczos bidiagonalization process with applications."
    SIAM Journal on Scientific Computing 21.6 (2000): 2257-2274.

[^Daniel1976]:
    Daniel, James W., et al. "Reorthogonalization and stable algorithms
    for updating the Gram-Schmidt QR factorization." Mathematics of
    Computation 30.136 (1976): 772-795.
```
"""
function extend!(
    log::ConvergenceHistory, A, L::PartialFactorization{T, Tr}, k::Int,
    orthleft::Bool=false, orthright::Bool=true, α::Real = 1/√2
    ) where {T,Tr}

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
    adjointA = adjoint(A)
    for j=l+1:k
        log.mtvps+=1
        mul!(q, adjointA, p) #q = A'p

        if orthright #Orthogonalize right Lanczos vector
            #Do double classical Gram-Schmidt reorthogonalization
            oldqnorm = norm(q)
            q -= L.Q*(L.Q'q)
            if norm(q) ≤ α * oldqnorm
                q -= L.Q*(L.Q'q)
            end
        end

        β = norm(q)
        q .*= inv(β)

        L.Q = [L.Q q]
        j==k && break

        log.mvps+=1
        #p = A*q - β*p
        mul!(p, A, q)
        p .-= β*view(L.P, :, j)

        if orthleft #Orthogonalize left Lanczos vector
            #Do double classical Gram-Schmidt reorthogonalization
            oldpnorm = norm(p)
            p -= L.P*(L.P'p)
            if norm(p) ≤ α * oldpnorm
                p -= L.P*(L.P'p)
            end
        end

        α = norm(p)
        p .*= inv(α)
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
    @assert Matrix(B) == [1 0 1; 0 2 2; 0 0 3]
end
