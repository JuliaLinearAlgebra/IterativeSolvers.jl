"""
    LimitedMemoryMatrix    

    A matrix which keeps only the last `memory` columns. Access outside these columns throws an error. Thus, the underlying storage is a matrix of size `(N, memory)` where `N` is the size of the first dimension
"""
mutable struct LimitedMemoryMatrix{T, V<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::V
    memory::Int
    hsize::Int
    
    function LimitedMemoryMatrix{T, V}(data::V, memory, hsize) where {T, V<:AbstractMatrix{T}}
        @assert Base.require_one_based_indexing(data)
        return new{T, V}(data, memory, hsize)
    end
end

function LimitedMemoryMatrix(mat::AbstractMatrix{T}, memory) where T
    data = similar(mat, size(mat, 1), memory)
    data[:, end-size(mat, 2)+1:end] .= mat
    return LimitedMemoryMatrix{T, typeof(data)}(data, memory, size(mat, 2))
end
function LimitedMemoryMatrix(vec::AbstractVector{T}, memory) where T
    data = similar(vec, size(vec, 1), memory)
    data[:, end] .= vec
    return LimitedMemoryMatrix{T, typeof(data)}(data, memory, 1)
end

Base.size(A::LimitedMemoryMatrix) = (size(A.data, 1), A.hsize)
Base.similar(A::LimitedMemoryMatrix) = LimitedMemoryMatrix(similar(A.data), A.memory, A.hsize)
Base.similar(A::LimitedMemoryMatrix, ::Type{T}) where T = LimitedMemoryMatrix(similar(A.data, T), A.memory, A.hsize)
Base.isempty(A::LimitedMemoryMatrix) = size(A, 1) == 0 || size(A, 2) == 0

Base.@propagate_inbounds function Base.setindex!(A::LimitedMemoryMatrix, v, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    lowest_stored_index = size(A, 2) - A.memory + 1
    if j < lowest_stored_index || j > size(A, 2)
        throw(
            ArgumentError(
                "Cannot set index ($i, $j) out of memory range, memory kept from $(max(1, lowest_stored_index)) to $(size(A, 2)))"
            )
        )
    else
        jj = j - lowest_stored_index + 1
        A.data[i, jj] = v
    end
    return v
end

Base.@propagate_inbounds function Base.getindex(A::LimitedMemoryMatrix, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    lowest_stored_index = size(A, 2) - A.memory + 1
    if j < lowest_stored_index || j > size(A, 2)
        throw(
            ArgumentError(
                "Cannot get index ($i, $j) out of memory range, memory kept from $(max(1, lowest_stored_index)) to $(size(A, 2)))"
            )
        )
    else
        @inbounds v = A.data[i, j - lowest_stored_index + 1]
    end
    return v
end

function Base.show(io::IO, ::MIME"text/plain", A::LimitedMemoryMatrix) where {T}
    print(io, Base.dims2string(size(A)), " ")
    Base.showarg(io, A.data, true)
end

function Base.show(io::IO, A::LimitedMemoryMatrix) where {T}
    print(io, Base.dims2string(size(A)), " ")
    Base.showarg(io, A.data, true)
end

"""
    hcat!(A::LimitedMemoryMatrix, B::AbstractVector)

Concatenates the vector `B` in-place to `A`, shifting the stored columns and increasing the
size of `A`
"""
function hcat!(A::LimitedMemoryMatrix, B::AbstractVector)
    @assert size(B, 1) == size(A, 1)
    A.hsize += 1
    A.data[:, 1:end-1] .= A.data[:, 2:end]
    A.data[:, end] .= B
    return A
end


"""
    LimitedMemoryUpperTriangularMatrix    

A matrix which keeps only the last `memory` columns. `setindex!` outside these columns
throws an error, otherwise entries are set to 0.
"""
mutable struct LimitedMemoryUpperTriangular{T, V<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::V
    memory::Int
    hsize::Int

    function LimitedMemoryUpperTriangular{T, V}(data::AbstractMatrix{T}, memory, hsize) where {T, V<:AbstractMatrix{T}}
        @assert Base.require_one_based_indexing(data)
        return new{T, typeof(data)}(data, memory, hsize)
    end
end

Base.size(A::LimitedMemoryUpperTriangular) = (A.hsize, A.hsize)
function LimitedMemoryUpperTriangular(mat::AbstractMatrix{T}, memory) where T
    data = similar(mat, memory, memory)
    fill!(data, 0)
    data[:, end-size(mat, 2)+1:end] .= mat
    return LimitedMemoryUpperTriangular{T, typeof(data)}(data, memory, size(mat, 2))
end
function LimitedMemoryUpperTriangular(vec::AbstractVector{T}, memory) where T
    data = similar(vec, memory, memory)
    fill!(data, 0)
    data[:, end] .= vec
    return LimitedMemoryUpperTriangular{T, typeof(data)}(data, memory, 1)
end
function LimitedMemoryUpperTriangular{T, V}(memory) where {T, V<:AbstractMatrix{T}}
    data = V(undef, memory, memory)
    fill!(data, 0)
    return LimitedMemoryUpperTriangular{T, typeof(data)}(data, memory, 0)
end

Base.@propagate_inbounds function Base.setindex!(A::LimitedMemoryUpperTriangular, v, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    lowest_stored_index = size(A, 2) - A.memory + 1
    if j < lowest_stored_index || (i > j) || (j-A.memory) >= i
        throw(
            ArgumentError(
                "Cannot set index ($i, $j) out of memory range, memory kept from $(max(1, lowest_stored_index)) to $(size(A, 2)))"
            )
        )
    else
        jj = j - lowest_stored_index + 1
        ii = i - lowest_stored_index + 1 - (size(A, 2) - j)
        A.data[ii, jj] = v
    end
    return v
end

Base.@propagate_inbounds function Base.getindex(A::LimitedMemoryUpperTriangular, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    lowest_stored_index_j = size(A, 2) - A.memory + 1
    if j < lowest_stored_index_j || (i > j) || (j-A.memory) >= i
        v = zero(eltype(A))
    else
        jj = j - lowest_stored_index_j + 1
        ii = i - lowest_stored_index_j + 1 + (size(A, 2) - j)
        @inbounds v = A.data[ii, jj]
    end
    return v
end

"""
    LimitedMemoryUpperHessenberg

A matrix which keeps only the last `memory` columns. `setindex!` outside these columns
throws an error, otherwise entries are set to 0.
"""
mutable struct LimitedMemoryUpperHessenberg{T, V<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::V
    memory::Int
    hsize::Int

    function LimitedMemoryUpperHessenberg{T, V}(data::AbstractMatrix{T}, memory, hsize) where {T, V<:AbstractMatrix{T}}
        @assert Base.require_one_based_indexing(data)
        return new{T, typeof(data)}(data, memory, hsize)
    end
end

Base.size(A::LimitedMemoryUpperHessenberg) = (A.hsize+1, A.hsize)
function LimitedMemoryUpperHessenberg(mat::AbstractMatrix{T}, memory) where T
    data = similar(mat, memory, memory)
    fill!(data, 0)
    data[end-memory+1:end, end-size(mat, 2)+1:end] .= mat
    return LimitedMemoryUpperHessenberg{T, typeof(data)}(data, memory, size(mat, 2))
end
function LimitedMemoryUpperHessenberg(vec::AbstractVector{T}, memory) where T
    data = similar(vec, memory, memory)
    fill!(data, 0)
    data[:, end] .= vec[end-memory+1:end]
    return LimitedMemoryUpperHessenberg{T, typeof(data)}(data, memory, 1)
end
function LimitedMemoryUpperHessenberg{T, V}(memory) where {T, V<:AbstractMatrix{T}}
    data = V(undef, memory, memory)
    fill!(data, 0)
    return LimitedMemoryUpperHessenberg{T, typeof(data)}(data, memory, 0)
end

Base.@propagate_inbounds function Base.setindex!(A::LimitedMemoryUpperHessenberg, v, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    lowest_stored_index = size(A, 2) - A.memory + 1
    if j < lowest_stored_index || (i > j+1) || (j+1-A.memory) >= i
        throw(
            ArgumentError(
                "Cannot set index ($i, $j) out of memory range, memory kept from $(max(1, lowest_stored_index)) to $(size(A, 2)))"
            )
        )
    else
        jj = j - lowest_stored_index + 1
        ii = i - lowest_stored_index - (size(A, 2) - j)
        A.data[ii, jj] = v
    end
    return v
end

Base.@propagate_inbounds function Base.getindex(A::LimitedMemoryUpperHessenberg, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)
    lowest_stored_index_j = size(A, 2) - A.memory + 1
    if j < lowest_stored_index_j || (i > j+1) || (j+1-A.memory) >= i
        v = zero(eltype(A))
    else
        jj = j - lowest_stored_index_j + 1
        ii = i - lowest_stored_index_j + (size(A, 2) - j)
        @inbounds v = A.data[ii, jj]
    end
    return v
end

function _grow_hcat!(A, v)
    if !isempty(A)
        @assert size(v, 1) == size(A, 1)+1
    end
    @assert sum(abs, v[1:end-A.memory]) < eps(real(eltype(A)))
    mem = min(A.memory, length(v))
    A.hsize += 1
    A.data[:, 1:end-1] .= A.data[:, 2:end]
    A.data[end-mem+1:end, end] .= v[end-mem+1:end]
    return A    
end