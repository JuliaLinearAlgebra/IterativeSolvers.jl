import Base: append!, eye, size

export ConvergenceHistory, KrylovSubspace

*(f::Function, v::AbstractVector)=f(v) #Syntax to mimic matvec product

type Eigenpair{S,T}
    val::S
    vec::Vector{T}
end

type ConvergenceHistory{T<:Real}
    isconverged::Bool
    threshold::T
    residuals::Vector{T}
end

type KrylovSubspace{T}
    A          #The linear operator that generates the subspace 
    n::Int     #Dimension of problem
    order::Int #Order (maximum size) of subspace
    v::Vector{Vector{T}} #The Krylov vectors
end

KrylovSubspace{T}(A, n::Int, order::Int=order, v::Vector{Vector{T}}=Vector{T}[])=
    KrylovSubspace(A, n, order, v)

KrylovSubspace{T}(A::AbstractMatrix{T}, order::Int=size(A,2), v::Vector{Vector{T}}=Vector{T}[])=
    KrylovSubspace(A, size(A,2), order, v)

#Reset an existing KrylovSubspace
function KrylovSubspace{T}(A::KrylovSubspace{T}, order::Int=size(A,2), v::Vector{Vector{T}}=Vector{T}[])
    A.order = order
    A.v = v
end

lastvec(K::KrylovSubspace) = K.v[end]
nextvec(K::KrylovSubspace) = K.A*lastvec(K)
size(K::KrylovSubspace) = length(K.v)
function size(K::KrylovSubspace, n::Int)
    if isa(K.A, AbstractMatrix)
        return size(K.A, n)
    else
        return n==2 ? size(K) : error("Illegal dimension $n of KrylovSubspace")
    end
end
eye{T}(K::KrylovSubspace{T}) = isa(K.A, Function) ? x->x : eye(T, size(K.A)...)

function append!{T}(K::KrylovSubspace{T}, w::Vector{T})
    if size(K) == K.order; shift!(K.v) end
    push!(K.v, w)
end
appendunit!{T}(K::KrylovSubspace{T}, w::Vector{T}) = append!(K, w/norm(w))

#Initialize the KrylovSubspace K with a random unit vector
function initrand!{T}(K::KrylovSubspace{T})
    v = convert(Vector{T}, randn(K.n))
    K.v = Vector{T}[v/norm(v)]
end

#Initialize the KrylovSubspace K with a specified nonunit vector
#If nonzero, try to normalize it
function init!{T}(K::KrylovSubspace{T}, v::Vector{T})
    K.v = Vector{T}[all(v.==zero(T)) ? v : v/norm(v)]
end

#Orthogonalize a vector v against the last p basis vectors defined by the
#KrylovSubspace K
function orthogonalize{T}(v::Vector{T}, K::KrylovSubspace{T}, p::Int;
    method::Symbol=:ModifiedGramSchmidt, normalize::Bool=false)
    Kk = K.v[max(1,end-p+1):end]
    if method == :GramSchmidt
        cs = T[dot(v, e) for e in Kk]
        for (i, e) in enumerate(Kk)
            v -= cs[i] * e
        end
    elseif method == :ModifiedGramSchmidt || method== :Householder
        #The numerical equivalence of ModifiedGramSchmidt and Householder was
        #established in doi:10.1137/0613015
        cs = zeros(T, p)
        for (i, e) in enumerate(Kk)
            cs[i] = dot(v, e)
            v-= cs[i] * e
        end
    else
        error("Unsupported orthogonalization method: $(method)")
    end
    normalize && (v /= convert(T, norm(v)))
    v, cs #Return orthogonalized vector and its coefficients
end

