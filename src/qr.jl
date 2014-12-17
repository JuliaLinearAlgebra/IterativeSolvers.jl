#Iterative algorithms for orthogonalization

import Base: Algorithm, qrfact!
export norm_none, norm_naive, norm_pythag,
       ReorthogonalizationAlg, never, always, rutishauser, giraudlangou,
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
    scale!(v, one(nrm)/nrm), nrm
end

#Pythagorean normalization
#Ref: doi:10.1007/s00211-006-0042-1
immutable norm_pythag <: NormalizationAlg end

function normalize!{T}(v::AbstractVector{T}, s::Real, p::Real, ::Type{norm_pythag})
    nrm = √(s-p)*√(s+p)
    scale!(v, one(nrm)/nrm), nrm
end

normalize!{T}(v::AbstractVector{T}, s::Real, u::AbstractVector{T}, ::Type{norm_pythag}) =
    normalize!(v, s, norm(u), norm_pythag)

const NORM_DEFAULT = norm_naive #Default is naive normalization

#######################
# Reorthogonalization #
#######################

#Define some commonly used criteria for checking whether to reorthogonalize

abstract Criterion
immutable never <: Criterion end
immutable always <: Criterion end

immutable rutishauser{T<:Real} <: Criterion
    K :: T
end
rutishauser()=rutishauser(10.0) #Rutishauser's original choice

immutable giraudlangou{T<:Real} <: Criterion
    K :: T
end

immutable ReorthogonalizationAlg{C<:Criterion} <: Algorithm
    criterion::C
    isforward::Bool
    maxiter::Int
end

ReorthogonalizationAlg(criterion::Criterion,
        isforward::Bool=true, maxiter::Integer=1) = 
    ReorthogonalizationAlg{C}(criterion, isforward, int(maxiter))
#Support calls like ReorthogonalizationAlg(never) and ReorthogonalizationAlg(always)
ReorthogonalizationAlg{C<:Criterion}(criterion::Type{C}) = 
    ReorthogonalizationAlg(criterion())
ReorthogonalizationAlg()=ReorthogonalizationAlg(never)

#Check if reorthogonalization is necessary and do it if it is
function reorthogonalize!(alg::ReorthogonalizationAlg, Q, R, k::Int)
    for iter = 1:alg.maxiter
        do_reorth(alg, Q, R, k) || continue
        for i = (alg.isforward ? (1:k-1) : (k-1:-1:1))
            d = Q[:,k]⋅Q[:,i]
            Q[:,k]-=d*Q[:,i]
            R[i,k]+=d
        end
    end
    nothing
end

function do_reorth(alg::ReorthogonalizationAlg, Q, R, k::Int)
    if isa(alg.criterion, never)
        return false
    elseif isa(alg.criterion, always)
        return true
    elseif isa(alg.criterion, rutishauser)
        #XXX This assumes that R[k,k] contains the original norm(Q[:,k])
        #    before orthogonalization!! Need to change rutishauser to
        #    stash this constant
        return R[k,k] > alg.criterion.K*norm(Q[:,k])
    elseif isa(alg.criterion, giraudlangou)
        return norm(R[1:k-1,k], 1) > alg.criterion.K*norm(Q[:,k])
    else
        throw(ArgumentError("Unknown ReorthogonalizationAlg: $alg"))
    end
end

################################
# Orthogonalization algorithms #
################################

#Type to store the output
if v"0.2" <= VERSION <= v"0.3-"
    immutable QRPair{T}
        Q::Matrix{T}
        R::Triangular{T}
        ϵ::Float64 #Residual
    end
else
    immutable QRPair{T}
        Q::Matrix{T}
        R::Triangular{T, Matrix{T}, :U}
        ϵ::Float64 #Residual
    end
end
QRPair(Q, R, ϵ::Real=0) = size(Q,2)==size(R,1)==size(R,2) ? QRPair(Q, R, float64(ϵ)) : throw(DimensionMismatch(""))

#The algorithm type hierarchy
abstract OrthogonalizationAlg <: Algorithm

#Classical Gram-Schmidt

immutable ClassicalGramSchmidtAlg{Norm<:NormalizationAlg, Reorth<:ReorthogonalizationAlg} <: OrthogonalizationAlg
    n::Union(Int, Nothing) #Number of Gram-Schmidt columns to compute
    normalize::Type{Norm}
    reorth::Reorth
end

typealias cgs ClassicalGramSchmidtAlg
cgs{Norm<:NormalizationAlg,Reorth<:ReorthogonalizationAlg}(n, norm::Type{Norm},
    reorth::Reorth=ReorthogonalizationAlg()) = cgs(n, norm, reorth)
cgs(n::Nothing) = cgs(n, NORM_DEFAULT)
cgs(n::Integer) = cgs(int(n), NORM_DEFAULT)
cgs() = cgs(nothing)

function qrfact!(Q, parameters::ClassicalGramSchmidtAlg)
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(eltype(Q), ncols, ncols)
    for k=1:ncols
        s=norm(Q[:,k])
        for i=1:k-1
            R[i,k] = Q[:,i]⋅Q[:,k]
        end
        for i=1:k-1
            Q[:,k]-= R[i,k]*Q[:,i]
        end
        reorthogonalize!(parameters.reorth, Q, R, k)
        Q[:,k], R[k,k] = normalize!(Q[:,k], parameters.normalize)
    end
    QRPair(Q, Triangular(R, :U))
end

function qrfact!(Q, parameters::ClassicalGramSchmidtAlg{norm_pythag})
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(eltype(Q), ncols, ncols)
    for k=1:ncols
        s=norm(Q[:,k])
        for i=1:k-1
            R[i,k] = Q[:,i]⋅Q[:,k]
        end
        for i=1:k-1
            Q[:,k]-= R[i,k]*Q[:,i]
        end
        reorthogonalize!(parameters.reorth, Q, R, k)
        Q[:,k], R[k,k] = normalize!(Q[:,k], s, R[1:k-1,k], norm_pythag)
    end
    QRPair(Q, Triangular(R, :U))
end

#Support qrfact!(A, cgs) call
qrfact!(A, ::Type{ClassicalGramSchmidtAlg}) = qrfact!(A, cgs()) 

#Columnwise modified Gram-Schmidt

immutable ColumnwiseModifiedGramSchmidtAlg{Norm<:NormalizationAlg, Reorth<:ReorthogonalizationAlg} <: OrthogonalizationAlg
    n::Union(Int, Nothing) #Number of Gram-Schmidt columns to compute
    normalize::Type{Norm}
    reorth::Reorth
end

typealias cmgs ColumnwiseModifiedGramSchmidtAlg
cmgs{Norm<:NormalizationAlg}(n, norm::Type{Norm},
    reorth::ReorthogonalizationAlg=ReorthogonalizationAlg()) = cmgs(n, norm, reorth)
cmgs(n::Nothing) = cmgs(n, NORM_DEFAULT)
cmgs(n::Integer) = cmgs(int(n), NORM_DEFAULT)
cmgs() = cmgs(nothing)

function qrfact!(Q, parameters::ColumnwiseModifiedGramSchmidtAlg)
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(eltype(Q), ncols, ncols)
    for k=1:ncols
    	for i=1:k-1
    	    R[i,k] = Q[:,i]⋅Q[:,k]
    	    Q[:,k]-= R[i,k]*Q[:,i]
    	end
        reorthogonalize!(parameters.reorth, Q, R, k)
        Q[:,k], R[k,k] = normalize!(Q[:,k], parameters.normalize)
    end
    QRPair(Q, Triangular(R, :U))
end

function qrfact!(Q, parameters::ColumnwiseModifiedGramSchmidtAlg{norm_pythag})
    ncols = parameters.n===nothing ? size(Q, 2) : parameters.n
    R=zeros(eltype(Q), ncols, ncols)
    for k=1:ncols
        s=norm(Q[:,k])
        for i=1:k-1
            R[i,k] = Q[:,i]⋅Q[:,k]
            Q[:,k]-= R[i,k]*Q[:,i]
        end
        reorthogonalize!(parameters.reorth, Q, R, k)
        Q[:,k], R[k,k] = normalize!(Q[:,k], s, R[1:k-1,k], norm_pythag)
    end
    QRPair(Q, Triangular(R, :U))
end

#Support qrfact!(A, cmgs) call
qrfact!(A, ::Type{ColumnwiseModifiedGramSchmidtAlg}) = qrfact!(A, cmgs()) 

include("../test/qr.jl")
