#Iterative algorithms for orthogonalization

import Base: Algorithm, qrfact!
export cgs, cmgs

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

immutable ClassicalGramSchmidtAlg <: OrthogonalizationAlg
    n::Union(Int, Nothing) #Number of Gram-Schmidt columns to compute
end

typealias cgs ClassicalGramSchmidtAlg
cgs() = cgs(nothing)

function qrfact!(Q, parameters::ClassicalGramSchmidtAlg)
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(ncols, ncols)
    for k=1:ncols
        for i=1:k-1
            R[i,k] = Q[:,i]⋅Q[:,k]
        end
        for i=1:k-1
            Q[:,k]-= R[i,k]*Q[:,i]
        end
        R[k,k] = norm(Q[:,k])
        Q[:,k]/= R[k,k]
    end
    QRPair(Q, Triangular(R, :U))
end
#Support qrfact!(A, cgs) call
qrfact!(A, ::Type{ClassicalGramSchmidtAlg}) = qrfact!(A, cgs()) 

#Columnwise modified Gram-Schmidt

immutable ColumnwiseModifiedGramSchmidtAlg <: OrthogonalizationAlg
    n::Union(Int, Nothing) #Number of Gram-Schmidt columns to compute
end

typealias cmgs ColumnwiseModifiedGramSchmidtAlg
cmgs() = cmgs(nothing)

function qrfact!(Q, parameters::ColumnwiseModifiedGramSchmidtAlg)
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(ncols, ncols)
    for k=1:ncols
	for i=1:k-1
	    R[i,k] = Q[:,i]⋅Q[:, k]
	    Q[:, k]-= R[i, k]*Q[:, i]
	end
	R[k,k] = norm(Q[:,k])
	Q[:,k]/= R[k,k]
    end
    QRPair(Q, Triangular(R, :U))
end
#Support qrfact!(A, cmgs) call
qrfact!(A, ::Type{ColumnwiseModifiedGramSchmidtAlg}) = qrfact!(A, cmgs()) 

include("../test/qr.jl")
