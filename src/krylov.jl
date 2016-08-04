import Base: append!, eye, size

export ConvergenceHistory, KrylovSubspace

"""
    KrylovSubspace

Collection of information on the Krylov subspace generated over the
course of an iterative Krylov solver.

**Fields**

* `A`: linear operator that generating the subspace.

* `n::Int`: dimension of problem.

* `order::Int`: order of subspace (maximum size).

* `v::Vector{Vector}`: Krylov vectors.

* `mvps::Int`: count of matrix-vector products.

**Constructors**

    KrylovSubspace(A, order)
    KrylovSubspace(A, n, order)
    KrylovSubspace(A, n, order, T)
    KrylovSubspace(A, n, order, v)

Create new `KrylovSubspace`.

#**Arguments**

* `T::Type`: type of the elements inside the Krylov vectors (eltype(v)).

* `v::Vector{Vector}`: initial Krylov subspace.


    KrylovSubspace(K, order=size(A,2), v=[])

Reset given `KrylovSubspace` `A`.

#**Arguments**

* `K::KrylovSubspace`: Krylov subspace to reset.

# Implements

* `Base`: `append!`, `eye`, `size`.

"""
type KrylovSubspace{T, OpT}
    A::OpT
    n::Int
    order::Int
    v::Vector{Vector{T}}
    mvps::Int
end
function KrylovSubspace(A, n::Int, order::Int, T::Type)
    KrylovSubspace{T,typeof(A)}(A, n, order, Array{T,1}[], 0)
end
function KrylovSubspace{T}(A, n::Int, order::Int, v::Vector{Vector{T}})
    KrylovSubspace{T,typeof(A)}(A, n, order, v, 0)
end
function KrylovSubspace(A, n::Int, order::Int)
    KrylovSubspace{eltype(A),typeof(A)}(A, n, order, Array{eltype(A),1}[], 0)
end
function KrylovSubspace(A, order::Int)
    KrylovSubspace{eltype(A),typeof(A)}(A, size(A,2), order, Array{eltype(A),1}[], 0)
end
function KrylovSubspace{T}(
        K::KrylovSubspace{T}, order::Int=size(K,2),
        v::Vector{Vector{T}}=Vector{T}[]
        )
    K.order = order
    K.v = v
    K.mvps = 0
end

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

"""
    lastvec(K)

Get last vector computed in the Krylov subspace `K`.

**Arguments**

* `K::KrylovSubspace`: Krylov subspace.

**Output**

* `::Vector`: last vector.

"""
lastvec(K::KrylovSubspace) = K.v[end]


"""
    nextvec(K)

Compute next vector in the Krylov subspace `K`.

**Arguments**

* `K::KrylovSubspace`: Krylov subspace.

**Output**

* `::Vector`: next vector.

"""
function nextvec(K::KrylovSubspace)
    K.mvps += 1
    K.A*lastvec(K)
end

"""
    appendunit!(K, w)

Append `normalize(w)` vector to krylob subspace `K`.

**Arguments**

* `K::KrylovSubspace`: Krylov subspace.

* `w::Vector`: vector to append.

"""
appendunit!{T}(K::KrylovSubspace{T}, w::Vector{T}) = append!(K, w/norm(w))

"""
    initrand!(K)

Initialize the Krylov subspace `K` with a random unit vector.

**Arguments**

* `K::KrylovSubspace`: Krylov subspace.

"""
function initrand!{T}(K::KrylovSubspace{T})
    K.v = Vector{T}[initrand!(Array(T, K.n))]
end

"""
    init!(K, v)

Initialize the KrylovSubspace K with a specified nonunit vector.

**Arguments**

* `K::KrylovSubspace`: Krylov subspace.

* `v::Vector`: initial vector.

"""
function init!{T}(K::KrylovSubspace{T}, v::Vector{T})
#    K.v = Vector{T}[all(v.==zero(T)) ? v : v/norm(v)]
    K.v = Vector{T}[copy(v)]
end

"""
    orthogonalize{T}(v, K, p)

Orthogonalize a vector `v` against the last `p` basis vectors defined by the
Krylov subspace `K`.

**Arguments**

* `v::Vector`: vector to orthogonalize.

* `K::KrylovSubspace`: subspace to orthogonalize `v` against.

* `p::Int=K.order`: last `p` remembered vectors.

#*Keywords*

* `method::Symbol=:ModifiedGramSchmidt`: orthogonalization method. Choose from:
    - `:GramSchmidt`: Gram Schmidt method.
    - `:ModifiedGramSchmidt`: Modified Gram Schmidt method.
    - `:Householder`: Householder method.

* `normalize::Bool=false`: normalize vector.

**Output**

* `::Vector`: orthogonalized vector.

* `::Vector`: orthogonalization coefficients.

"""
function orthogonalize{T}(v::Vector{T}, K::KrylovSubspace{T}, p::Int=K.order;
    method::Symbol=:ModifiedGramSchmidt, normalize::Bool=false)
    Kk = K.v[max(1,end-p+1):end]
    if method == :GramSchmidt
        cs = T[dot(e, v) for e in Kk]
        for (i, e) in enumerate(Kk)
            v -= cs[i] * e
        end
    elseif method == :ModifiedGramSchmidt
        cs = zeros(T, p)
        for (i, e) in enumerate(Kk)
            cs[i] = dot(e, v)
            v-= cs[i] * e
        end
    elseif method == :Householder
        #Use brute force Householder
        Q, R = qr([Kk v])
        cs = R[:, end]
        v = Q[:, end]
    else
        error("Unsupported orthogonalization method: $(method)")
    end
    normalize && (v /= convert(T, norm(v)))
    v, cs
end

"""
    init!(a, K)

Compute a linear combination of Krylov vectors.

**Arguments**

* `a::Vector`: initial vector.

* `K::KrylovSubspace`: Krylov subspace.

**Output**

* `::Vector`: vector in the Krylov basis.

"""
function *(a::AbstractVector, K::KrylovSubspace)
    N = length(a)
    N==length(K.v) && throw(DimensionMismatch("Tried to compute linear combination of $N Krylov vectors but only $(length(K.v)) are present."))
    w = a[N] * K.v[N]
    for j = (N-1):-1:1
        #Index 1 = oldest remembered vector
        w += a[j] * K.v[j]
    end
    w
end
