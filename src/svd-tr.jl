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


function thickrestartbidiag(A, q::AbstractVector, k::Int, l::Int;
    maxiter::Int=10, tol::Real=√eps())

    @assert k>l
    L = build(A, q, k)

    local F
    for i=1:maxiter
        info("Iteration $i")
        #@assert size(L.B) == (k, k)
        F = svdfact(L.B)
        L = truncate!(A, L, F, l)
        extend!(A, L, k)
        isconverged(L, F, tol) && break
    end
    F[:S][1:l], L
end

function isconverged(L::BidiagonalFactorization,
        F::Base.LinAlg.SVD, tol::Real)

    σ = F[:S]
    Δσ= L.β*abs(F[:U][end, :])

    for i in eachindex(σ)
        if Δσ[i]<√tol
            println(σ[i], " ± ", Δσ[i], " ")
        end
    end
    all(Δσ .< tol)
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
    p = Array(T, m)
    for j=1:k
        #p = A*q
        A_mul_B!(p, A, q)
        if j>1
            βs[j-1] = β
            Base.LinAlg.axpy!(-β, sub(P, :, j-1), p) #p -= β*P[:,j-1]
        end
        α = norm(p)
        p[:] = p/α

        αs[j] = α
        P[:, j] = p

        Ac_mul_B!(q, A, p) #q = A'p
        Base.LinAlg.axpy!(-1.0, sub(Q, :, 1:j)*(sub(Q, :, 1:j)'q), q) #orthogonalize

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
    @assert size(L.P) == (m, k)
    @assert size(L.Q) == (n, k+1)

    @assert all([norm(L.Q[:,i]) ≈ 1 for i=1:size(L.Q,2)])
    L.Q = [sub(L.Q, :,1:k)*sub(F[:V], :,1:l) sub(L.Q, :, k+1)]

    #Be pedantic about ensuring normalization
    #L.Q = qr(L.Q)[1]
    @assert all([norm(L.Q[:,i]) ≈ 1 for i=1:size(L.Q,2)])

    f = A*sub(L.Q, :, l+1)
    ρ = L.β * reshape(F[:U][end, 1:l], l)
    L.P = sub(L.P, :, 1:k)*sub(F[:U], :, 1:l)

    #@assert ρ[i] ≈ f⋅L.P[:, i]
    f -= L.P*ρ
    α = norm(f)
    f[:] = f/α
    L.P = [L.P f]

    g = A'f - α*L.Q[:, end]
    β = norm(g)
    g[:] = g/β
    L.β = β
    L.B = BrokenArrowBidiagonal([F[:S][1:l]; α], ρ, typeof(β)[])

    L
end

function extend!(A, L::BidiagonalFactorization, k::Int)
    l = size(L.P, 2)-1
    p = L.P[:, end]

    local β
    m, n = size(A)
    q = zeros(n)
    for j=l+1:k
        Ac_mul_B!(q, A, p) #q = A'p
        q -= L.Q*(L.Q'q)   #orthogonalize
        β = norm(q)
        q[:] = q/β
        push!(L.B.ev, β)

        L.Q = [L.Q q]
        j==k && break

        #p = A*q - β*p
        A_mul_B!(p, A, q)
        Base.LinAlg.axpy!(-β, sub(L.P, :, j), p)

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

let
    m = 300
    n = 200
    k = 10
    l = 5

    A = randn(m,n)
    q = randn(n)|>x->x/norm(x)
    σ, L = thickrestartbidiag(A, q, k, l)

    @assert norm(σ - svdvals(A)[1:l]) < 0.001
end

using Base.Profile
let
    srand(1)
    m = 30000
    n = 20000
    k = 200
    l = 100

    A = sprandn(m,n,0.1)
    q = randn(n)|>x->x/norm(x)
    Profile.clear()
    @time @profile L = thickrestartbidiag(A, q, k, l; tol=1e-4)
    Profile.print()
end
