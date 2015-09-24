import Base.LinAlg: axpy!

export svdvals_tr

"""
Matrix of the form
 d1          a_1
    d2       a_2
      ...    ...
         d_l a_l
             d_l+1 e_l+1
                   ...  e_k-1
                        d_k
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
Partial factorization object which is an approximation of a matrix

    A ≈ P * [B 0; 0 β] * Q
"""
type PartialFactorization{T, Tr} <: Factorization{T}
    P :: Matrix{T}
    Q :: Matrix{T}
    B :: AbstractMatrix{Tr}
    β :: Tr
end


"""
The thick-restarted variant of Golub-Kahan-Lanczos bidiagonalization

# Inputs

- `A` : The matrix or matrixlike object whose singular values are desired
- `l` : The number of singular values requested.
        Default: 6

# Keyword inputs

- `v0` : The starting guess vector in the range of `A`.
         The length of `q` should be the number of columns in `A`.
         Default: A random unit vector.
- `k` : The maximum number of Lanczos vectors to compute before restarting.
        Default: `2*l`
- `j` : The number of vectors to keep at the end of the restart.
        Default: `l`. We don't recommend j < l.
- `maxiter`: Maximum number of iterations to run
             Default: `minimum(size(A))`
- `verbose`: Whether to print information at each iteration
             Default: false
- `tol`    : Maximum absolute error in each desired singular value.
             Default: `√eps()`
- `reltol` : Maximum error in each desired singular value relative to the
             estimated norm of the input matrix.  Default: `√eps()`
- `restart`: Which restarting algorithm to use. Valid choices are:
             - `:ritz`: Thick restart with Ritz values [Wu2000] (default)
             - `:harmonic`: Restart with harmonic Ritz values [Baglama2005]

- `doplot` : Plot a history of the Ritz value convergence. Requires the
             UnicodePlots.jl package to be installed. Default: false

# Output

- `Σ`: A list of the desired singular values
- `L`: The computed partial factorization of A


# Implementation notes

The implementation of thick restarting follows closely that of SLEPc as
described in [Hernandez2008].

# References

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

"""
function svdvals_tr(A, l::Int=6; k::Int=2l,
    j::Int=l, v0::AbstractVector = Vector{eltype(A)}(randn(size(A, 2))) |> x->scale!(x, inv(norm(x))),
    maxiter::Int=minimum(size(A)), tol::Real=√eps(), reltol::Real=√eps(),
    verbose::Bool=false, method::Symbol=:ritz, doplot::Bool=false, vecs=:none)

    T0 = time_ns()
    @assert k>l
    L = build(A, v0, k)

    #Save history of Ritz values
    ritzvalhist = []
    convhist = []

    if doplot
       ttyh, ttyw = Base.tty_size()
       ttyh -= l+1+5+2
       ttyw -= 12
    end

    local F
    for iter in 1:maxiter
        #@assert size(L.B) == (k, k)
        F = svdfact(L.B)
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
        conv = isconverged(L, F, l, tol, reltol, verbose)

        push!(convhist, conv)
        push!(ritzvalhist, F[:S])

        if doplot
            if isdefined(Main, :UnicodePlots)
                xs = convert(Vector{Vector{Int}}, map(x->fill(x[1], length(x[2])), enumerate(ritzvalhist)))
                layer1 = Main.UnicodePlots.scatterplot([xs...;], [ritzvalhist...;], height=ttyh, width=ttyw)
                display(Main.UnicodePlots.scatterplot!(layer1,
                    [[fill(x[1], sum(conv)) for x in enumerate(convhist)]...;],
                    [[x[1:l][conv] for x in ritzvalhist]...;],
                    height=ttyh, width=ttyw, color=:red))
            else
                warn("UnicodePlots not found; no plotsies")
            end
        end

        #Lock
        if method == :ritz
            for i in eachindex(conv)
                if conv[i]
                    L.B.av[i] = 0
                end
            end
        end
        all(conv) && break
    end
    F[:S][1:l], L
end

"""
# References

The simple error bound dates back at least to Wilkinson's classic book
[Wilkinson1965:Ch.3 §53 p.170]

@book{Wilkinson1965,
    address = {Oxford, UK},
    author = {J H Wilkinson},
    publisher = {Oxford},
    title = {The Algebraic Eigenvalue Problem},
    year = 1965
}

The Rayleigh-Ritz bounds were presented in [Wilkinson1965:Ch.3 §54-55 p.173,
    Yamamoto1980, Ortega1990]

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

"""
function isconverged(L::PartialFactorization,
        F::Base.LinAlg.SVD, k::Int, tol::Real, reltol::Real, verbose::Bool=false)

    @assert tol ≥ 0

    σ = F[:S][1:k]
    Δσ= L.β*abs(F[:U][end, 1:k])

    #Best available eigenvalue bounds
    δσ = copy(Δσ)

    #Update better error bounds from Rayleigh-Ritz
    #Reference: Ortega, 1972
    #In theory these only apply if we know all the values
    #But so long as the gap between the converged Ritz values and all the
    #others is larger than d then we're fine.
    #Can also get some eigenvector statistics!!
    let
        d = Inf
        for i in eachindex(σ), j=1:i-1
            d = min(d, abs(σ[i] - σ[j]))
        end
        verbose && println("Smallest empirical spectral gap: ", d)
        #Chatelein 1993 - normwise backward error associated with approximate invariant subspace
        #Use the largest singular (Ritz) value to estimate the 2-norm of the matrix
        #verbose && println("Normwise backward error associated with subspace: ", L.β/σ[1])
        for i in eachindex(Δσ)
            α = Δσ[i]

            #Simple error bound
            #if 2α ≤ d

            if 2α ≤ d
                #Rayleigh-Ritz bounds
                x = α/(d-α)*√(1+(α/(d-α))^2)
                x=abs(x)
                #verbose && println("Rayleigh-Ritz error bound on eigenvector: $x")
                #2α ≤ d && (δσ[i] = min(δσ[i], x))

                y = α^2/d #[Wilkinson:Ch.3 Appendix (4), p.188]
                #verbose && println("Rayleigh-Ritz error bound on eigenvalue: $y")
                2α ≤ d && (δσ[i] = min(δσ[i], y))
            end

	    verbose && println("Ritz value ", i, ": ", σ[i], " ± ", signif(δσ[i], 3))

            #Estimate of the normwise backward error [Deif 1989]
            ##Use the largest singular (Ritz) value to estimate the 2-norm of the matrix
            #verbose && println("Normwise backward error estimate: ", α/σ[1])

            #Geurts, 1982 - componentwise backward error also known
            #A. J. Geurts, (1982), A contribution to the theory of condition,
            #Numer.  Math., #39, 85-96.
        end
    end

    #Estimate condition number and see if two-sided reorthog is needed
    if verbose && (F[:S][1]/F[:S][end]) > 1/√eps()
        warn("Two-sided reorthogonalization should be used but is not implemented")
    end

    conv = δσ[1:k] .< max(tol, reltol*σ[1])
end

#Hernandez2008
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
Thick restart (with ordinary Ritz values)

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
Thick restart with harmonic Ritz values

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

    F2 = svdfact!([diagm(F0[:S]) ρ'], thin=false)

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

    R = sub(R,1:k+1,1:k) + sub(R,:,k+1)*Mend

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

#Hernandez2008
function extend!{T,Tr}(A, L::PartialFactorization{T, Tr}, k::Int)
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
        q -= L.Q*(L.Q'q)   #orthogonalize
        β = norm(q)
        scale!(q, inv(β))

        L.Q = [L.Q q]
        j==k && break

        #p = A*q - β*p
        A_mul_B!(p, A, q)
        Base.LinAlg.axpy!(-β, sub(L.P, :, j), p)

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


