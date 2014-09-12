import Base: start, next, done

#The usual Krylov space machinery
type Krylov{T}
    A
    v₀ :: AbstractVector{T}
end

#Termination criteria
type Terminator{T}
    tol :: T
    maxiter :: Int
end

# __Definition 2.1.__ The Arnoldi factorization is a truncated Hessenberg factorization of a matrix $A$ such that
#
# $$ AV = VH + r e_k^T $$
#
# where $V \in \mathbb{R}^{n\times k}$, $V^T V = I_k$, $H \in \mathbb{R}^{k\times k}$ is upper Hessenberg, $r \in \mathbb{R}^n$ with $0=V^T r$.

type ArnoldiFact{T} <: Factorization{T}
    V :: Matrix{T}
    H :: Matrix{T}
    r :: Vector{T}
end

function ArnoldiFact{T}(V::Matrix{T}, H::Matrix{T}, r::Vector{T})
    #Check dimensions
    @assert size(V,2) == size(H,1)
    @assert size(V,1) == size(r,1)
    ArnoldiFact{T}(V, H, r)
end

#The Arnoldi iterator

abstract Factorizer

type Arnoldi{T} <: Factorizer
    K :: Krylov{T}
    term :: Terminator{T}
end

function start{T}(P::Arnoldi{T})
    r = P.K.v₀
    n = size(r, 1)
    V = zeros(T,n,0)
    H = zeros(T,0,0)
    ArnoldiFact{T}(V, H, r)
end

# To write out the algorithms later, it's also convenient to define this `UnitVec` type which represents the canonical basis vector $e_k$. The nice thing is that pre-/post-multiplying $A$ by $e_k$ is equivalent to extracting the $k$th row or column of $A$ respectively.
type UnitVec
    k :: Int
end

*(A::AbstractArray, e::UnitVec) = A[:,e.k]
*(e::UnitVec, A::AbstractArray) = A[e.k,:]

# __The `Arnoldi` function (Algorithm 3.7; Sorensen, 1992)__.
#
# Input: $AV - VH = re_k^T$ with $V^T V = I_k, V^T r = 0$
#
# Output: $AV - VH = re_{k+p}^T$ with $V^T V = I_{k+p}, V^T r = 0$

function next{T}(P::Arnoldi{T}, F::ArnoldiFact{T})
    V, H, r = F.V, F.H, F.r
    β = norm(r)

    H = size(H,2)==0 ? zeros(1,0) :
      [H; [zeros(1,size(H,2)-1) β]]
    v = r / β
    V = [V v]

    w = K.A*v

    h = V'w
    H = [H h]

    r = w - V*h
    #Gram-Schmidt
    for i=1:2
      s = V'r
      r-= V*s
      h+= s
      norm(s) < P.term.tol*norm(r) && break
    end
    @assert abs(norm(K.A*V-V*H) - norm(r)) < P.term.tol
    F  = ArnoldiFact(V, H, r)
    F, F
end

function done{T}(P::Arnoldi{T}, F::ArnoldiFact{T})
    #TODO maxiter
    β = norm(F.r)
    β < P.term.tol
end

#Test of Arnoldi iterator

n=10
M=randn(n,n)#;M+=M'
K = Krylov(M, randn(n))
for (iter, Ar) in enumerate(Arnoldi(K, Terminator(1e-9, n)))
    println("Iteration $iter:\t residual norm = ", norm(Ar.r))
    @assert abs(norm(K.A*Ar.V-Ar.V*Ar.H) - norm(Ar.r)) < sqrt(eps())
end


# # References
#
# 1. D. C. Sorensen, "Implicit application of polynomial filters in a $k$-step Arnoldi method", _SIAM J. Matrix Anal. Appl._ 13 (1), pp. 357-385, January 1992. [doi:10.1137/0613025](http://epubs.siam.org/doi/abs/10.1137/0613025)
