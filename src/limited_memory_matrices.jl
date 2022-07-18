"""
    LimitedMemoryMatrix    

    A matrix which keeps only the last `memory` columns. Access outside these columns throws an error. Thus, the underlying storage is a matrix of size `(N, memory)` where `N` is the size of the first dimension
"""
mutable struct LimitedMemoryMatrix{T, V<:AbstractMatrix{T}} <: AbstractMatrix{T}
    value::V
    memory::Int
    hsize::Int
    
    function LimitedMemoryMatrix{T, V}(value::V, memory, hsize) where {T, V<:AbstractMatrix{T}}
        @assert Base.require_one_based_indexing(value)
        return new{T, V}(value, memory, hsize)
    end
end

function LimitedMemoryMatrix(mat::AbstractMatrix{T}, memory) where T
    value = similar(mat, size(mat, 1), memory)
    value[:, end-size(mat, 2)+1:end] .= mat
    return LimitedMemoryMatrix{T, typeof(value)}(value, memory, size(mat, 2))
end
function LimitedMemoryMatrix(vec::AbstractVector{T}, memory) where T
    value = similar(vec, size(vec, 1), memory)
    value[:, end] .= vec
    return LimitedMemoryMatrix{T, typeof(value)}(value, memory, 1)
end

Base.size(A::LimitedMemoryMatrix) = (size(A.value, 1), A.hsize)
Base.similar(A::LimitedMemoryMatrix) = LimitedMemoryMatrix(similar(A.value), A.memory, A.hsize)
Base.similar(A::LimitedMemoryMatrix, ::Type{T}) where T = LimitedMemoryMatrix(similar(A.value, T), A.memory, A.hsize)
Base.isempty(A::LimitedMemoryMatrix) = size(A, 1) == 0 || size(A, 2) == 0

Base.@propagate_inbounds function Base.setindex!(A::LimitedMemoryMatrix, v, i::Integer, j::Integer)
    lowest_stored_index = size(A, 2) - A.memory + 1
    if j < lowest_stored_index || j > size(A, 2)
        throw(
            ArgumentError(
                "Cannot set index ($i, $j) out of memory range, memory kept from $(max(1, lowest_stored_index)) to $(size(A, 2)))"
            )
        )
    else
        jj = j - lowest_stored_index + 1
        A.value[i, jj] = v
    end
    return v
end

Base.@propagate_inbounds function Base.getindex(A::LimitedMemoryMatrix, i::Integer, j::Integer)
    lowest_stored_index = size(A, 2) - A.memory + 1
    if j < lowest_stored_index || j > size(A, 2)
        throw(
            ArgumentError(
                "Cannot get index ($i, $j) out of memory range, memory kept from $(max(1, lowest_stored_index)) to $(size(A, 2)))"
            )
        )
    else
        @inbounds v = A.value[i, j - lowest_stored_index + 1]
    end
    return v
end

function Base.show(io::IO, ::MIME"text/plain", A::LimitedMemoryMatrix) where {T}
    print(io, Base.dims2string(size(A)), " ")
    Base.showarg(io, A.value, true)
end

"""
    hcat!(A::LimitedMemoryMatrix, B::AbstractVector)

Concatenates the vector `B` in-place to `A`, shifting the stored columns and increasing the
size of `A`
"""
function hcat!(A::LimitedMemoryMatrix, B::AbstractVector)
    @assert size(B, 1) == size(A, 1)
    A.hsize += 1
    A.value[:, 1:end-1] .= A.value[:, 2:end]
    A.value[:, end] .= B
    return A
end

"""
    LimitedMemoryBandedMatrix    

A matrix which keeps only the last `memory` columns. Access outside these columns throws an
  error. Additionally, the matrix is understood to be banded, such that only values
  `band_offset` from the diagonal are non-zero, where the band is of width `band_width`. Thus, the underlying storage is a matrix of size `(band_width, memory)`
"""
mutable struct LimitedMemoryBandedMatrix{T, V<:AbstractMatrix{T}} <: AbstractMatrix{T}
    value::V
    memory::Int
    size::Tuple{Int, Int}
    band_offset::Int
    band_width::Int
end