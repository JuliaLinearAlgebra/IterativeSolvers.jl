import Base.size
import Base.LinAlg.BlasFloat
import Base.append!

type KrylovSubspace{T}
    A::AbstractMatrix{T} #The matrix that generates the subspace 
    n::Int               #Dimension of problem
    maxdim::Int          #Maximum size of subspace
    v::Vector{Vector{T}} #The Krylov vectors
end
function KrylovSubspace{T}(A::AbstractMatrix{T}, maxdim::Int)
    n = size(A, 2)
    v = Vector{T}[]
    K = KrylovSubspace(A, n, maxdim, v)
end
KrylovSubspace(A::AbstractMatrix) = KrylovSubspace(A, size(A, 2))

lastvec(K::KrylovSubspace) = K.v[end]
nextvec(K::KrylovSubspace) = K.A * K.v[end]
size(K::KrylovSubspace) = length(K.v)

function append!{T}(K::KrylovSubspace{T}, w::Vector{T})
    if size(K) == K.maxdim; shift!(K.v) end
    push!(K.v, w)
end
appendunit!{T}(K::KrylovSubspace{T}, w::Vector{T}) = append!(K, w/norm(w))

#Initialize the KrylovSubspace K with a random unit vector
function initrand!{T<:Real}(K::KrylovSubspace{T})
    v = convert(Vector{T}, randn(K.n))
    K.v = Vector{T}[v/norm(v)]
end

#Orthogonalize a vector v against the basis defined by the KrylovSubspace K
function orthogonalize!{T}(v::Vector{T}, K::KrylovSubspace{T}; method::Symbol=:GramSchmidt)
    if method == :GramSchmidt
        for e in K.v
            v -= dot(v, e) * e
        end
    else
        error(string("Unsupported orthgonalization method: ", method))
    end
    v /= norm(v) #Always return unit vector
end

