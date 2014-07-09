#Iterative algorithms for orthogonalization

import Base: Algorithm, qrfact!
export norm_none, norm_naive, norm_pythag,
       cgs, cmgs

#########################
# Normalization methods #
#########################

abstract NormalizationAlg <: Algorithm
normalize!(v) = normalize!(v, NORM_DEFAULT)

#No normalization
immutable norm_none <: NormalizationAlg end
normalize!(v::AbstractVector, ::Type{norm_none}) = v, 1.0

#Naive normalization
immutable norm_naive <: NormalizationAlg end
function normalize!(v::AbstractVector, ::Type{norm_naive})
    nrm = norm(v)
    v/nrm, nrm
end

#Pythagorean normalization
#Ref: doi:10.1007/s00211-006-0042-1
immutable norm_pythag <: NormalizationAlg end

function normalize!{T}(v::AbstractVector{T}, s::Real, p::Real, ::Type{norm_pythag})
    nrm = √(s-p)*√(s+p)
    v/nrm, nrm
end

normalize!{T}(v::AbstractVector{T}, s::Real, u::AbstractVector{T}, ::Type{norm_pythag}) =
    normalize!(v, s, norm(u), norm_pythag)

const NORM_DEFAULT = norm_naive #Default is naive normalization

################################
# Orthogonalization algorithms #
################################

#Type to store the output
immutable QRPair{T}
    Q::Matrix{T}
    R::Triangular{T, Matrix{T}, :U}
    ϵ::Float64 #Residual
end
QRPair(Q, R, ϵ::Real=0) = size(Q,2)==size(R,1)==size(R,2) ? QRPair(Q, R, float64(ϵ)) : throw(DimensionMismatch(""))

#The algorithm type hierarchy
abstract OrthogonalizationAlg <: Algorithm

#Classical Gram-Schmidt

immutable ClassicalGramSchmidtAlg{A<:NormalizationAlg} <: OrthogonalizationAlg
    n::Union(Int, Nothing) #Number of Gram-Schmidt columns to compute
    normalize::Type{A}
end

typealias cgs ClassicalGramSchmidtAlg
cgs(n::Nothing) = cgs(n, NORM_DEFAULT)
cgs(n::Integer) = cgs(int(n), NORM_DEFAULT)
cgs() = cgs(nothing)

function qrfact!(Q, parameters::ClassicalGramSchmidtAlg)
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(ncols, ncols)
    for k=1:ncols
        parameters.normalize===norm_pythag && (s=norm(Q[:,k]))
        for i=1:k-1
            R[i,k] = Q[:,i]⋅Q[:,k]
        end
        for i=1:k-1
            Q[:,k]-= R[i,k]*Q[:,i]
        end
        Q[:,k], R[k,k] = parameters.normalize===norm_pythag ?
            normalize!(Q[:,k], s, R[1:k-1,k], norm_pythag) :
            normalize!(Q[:,k], parameters.normalize)
    end
    QRPair(Q, Triangular(R, :U))
end
#Support qrfact!(A, cgs) call
qrfact!(A, ::Type{ClassicalGramSchmidtAlg}) = qrfact!(A, cgs()) 

#Columnwise modified Gram-Schmidt

immutable ColumnwiseModifiedGramSchmidtAlg{A<:NormalizationAlg} <: OrthogonalizationAlg
    n::Union(Int, Nothing) #Number of Gram-Schmidt columns to compute
    normalize::Type{A}
end

typealias cmgs ColumnwiseModifiedGramSchmidtAlg
cmgs(n::Nothing) = cmgs(n, NORM_DEFAULT)
cmgs(n::Integer) = cmgs(int(n), NORM_DEFAULT)
cmgs() = cmgs(nothing)

function qrfact!(Q, parameters::ColumnwiseModifiedGramSchmidtAlg)
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(ncols, ncols)
    for k=1:ncols
        parameters.normalize===norm_pythag && (s=norm(Q[:,k]))
    	for i=1:k-1
    	    R[i,k] = Q[:,i]⋅Q[:, k]
    	    Q[:, k]-= R[i, k]*Q[:, i]
    	end
        Q[:,k], R[k,k] = parameters.normalize===norm_pythag ?
            normalize!(Q[:,k], s, R[1:k-1,k], norm_pythag) :
            normalize!(Q[:,k], parameters.normalize)
    end
    QRPair(Q, Triangular(R, :U))
end
#Support qrfact!(A, cmgs) call
qrfact!(A, ::Type{ColumnwiseModifiedGramSchmidtAlg}) = qrfact!(A, cmgs()) 

include("../test/qr.jl")
