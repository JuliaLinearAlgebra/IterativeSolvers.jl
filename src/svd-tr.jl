#Matrix of the form
# d1          a_1
#    d2       a_2
#      ...    ...
#         d_l a_l
#             d_l+1 e_l+1
#                   ...  e_k-1
#                        d_k
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

Base.svdfact(B::BrokenArrowBidiagonal) = svdfact(full(B))

type BidiagonalFactorization{T, Tr} <: Factorization{T}
    P :: T
    Q :: T
    B :: Union(Bidiagonal{Tr}, BrokenArrowBidiagonal{Tr})
    β :: Tr
end


function thickrestartbidiag(A, q::AbstractVector, k::Int, l::Int, maxiter::Int=10)
    @assert k>l
    L = build(A, q, k)
    for i=1:maxiter
        info("Iteration $i")
        @show size(L.B), k
        @assert size(L.B) == (k, k)
        F = svdfact(L.B)
        @assert all([norm(F[:V][:,i]) ≈ 1 for i in 1:size(F[:V], 2)])
        @assert all([norm(F[:U][:,i]) ≈ 1 for i in 1:size(F[:U], 2)])
        L = truncate!(A, L, F, l)
        extend!(A, L, k)
        isconverged(L, F) && break
    end
    L
end

function isconverged(L::BidiagonalFactorization,
        F::Base.LinAlg.SVD)

    σ = F[:S]
    Δσ= L.β*abs(F[:U][end, :])

    for i in eachindex(σ)
        if Δσ[i]<√1e-7
            println(σ[i], " ± ", Δσ[i], " ")
        end
    end
    all(Δσ .< 1e-7)
end

function build{T}(A, q::AbstractVector{T}, k::Int)
    m, n = size(A)
    Tr = typeof(real(one(T)))

    P = Array(T, m, k)
    Q = Array(T, n, k+1)
    αs= Array(Tr, k)
    βs= Array(Tr, k-1)
    Q[:, 1] = q
    local β
    for j=1:k
        p = A*q
        if j>1
            βs[j-1] = β
            p -= β*P[:,j-1]
        end
        α = norm(p)
        p[:] = p/α

        αs[j] = α
        P[:, j] = p

        q = A'p
        for i=1:j
            γ = Q[:,i]⋅q
            q-= γ*Q[:,i]
        end
        β = norm(q)
        q[:] = q/β
        Q[:,j+1] = q
    end
    BidiagonalFactorization(P, Q, Bidiagonal(αs, βs, true), β)
end

#By default, truncate to l largest singular triplets
function truncate!(A, L::BidiagonalFactorization,
        F::Base.LinAlg.SVD, l::Int)

    k = size(F[:V], 1)
    m, n = size(A)
    @show size(L.P), (m, k)
    @assert size(L.P) == (m, k)
    @assert size(L.Q) == (n, k+1)

    @assert all([norm(L.Q[:,i]) ≈ 1 for i=1:size(L.Q,2)])
    L.Q = [L.Q[:,1:k]*F[:V][:,1:l] L.Q[:, end]]

    #Be pedantic about ensuring normalization
    #L.Q = qr(L.Q)[1]
    @assert all([norm(L.Q[:,i]) ≈ 1 for i=1:size(L.Q,2)])

    f = A*L.Q[:, end]
    ρ = L.β * reshape(F[:U][end, 1:l], l)
    L.P = L.P[:, 1:k]*F[:U][:,1:l]
    for i=1:l
        #if !(ρ[i] ≈ f⋅L.P[:, i])
        #    @show ρ[i], f⋅L.P[:, i]
        #end
        #ρ[i] = f⋅L.P[:, i]
        f -= ρ[i]*L.P[:, i]
    end
    α = norm(f)
    f[:] = f/α
    L.P = [L.P f]

    g = A'f - α*L.Q[:, end]
    β = norm(g)
    g[:] = g/β
    #L.Q = [L.Q g]
    L.β = β
    L.B = BrokenArrowBidiagonal([F[:S][1:l]; α], ρ, typeof(β)[])

    L
end

function extend!(A, L::BidiagonalFactorization, k::Int)
    l = size(L.P, 2)-1
    p = L.P[:, end]
    #for i=1:l
    #    @show i, norm(p)
    #    @assert norm(L.P[:,i]) ≈ 1
    #    ρ = p⋅L.P[:, i]
    #    p -= ρ*L.P[:, i]
    #end

    local β
    for j=l+1:k
        q = A'p
        for i=1:j
            @assert norm(L.Q[:,i]) ≈ 1
            γ = L.Q[:,i]⋅q
            q -= γ*L.Q[:,i]
        end
        β = norm(q)
        q[:] = q/β
        push!(L.B.ev, β)

        L.Q = [L.Q q]
        j==k && break
        p = A*q - β*p
        α = norm(p)
        p[:] = p/α
        push!(L.B.dv, α)
        L.P = [L.P p]
    end
    L.β = β
    L
end

let
    B = BrokenArrowBidiagonal([1, 2, 3], [1, 2], Int[])
    @assert full(B) == [1 0 1; 0 2 2; 0 0 3]
end

#A simple test
let
    m = 3000
    n = 2000
    k = 20
    l = 10

    A = randn(m,n)
    q = randn(n)|>x->x/norm(x)
    @time L = thickrestartbidiag(A, q, k, l)

    @show svdvals(A)[1:l]
end
