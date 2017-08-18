function laplace_matrix(::Type{T}, n, dims) where T
    D = second_order_central_diff(T, n);
    A = copy(D);

    for idx = 2 : dims
        A = kron(A, speye(n)) + kron(speye(size(A, 1)), D);
    end

    A
end

second_order_central_diff(::Type{T}, dim) where {T} = convert(SparseMatrixCSC{T, Int}, SymTridiagonal(fill(2 * one(T), dim), fill(-one(T), dim - 1)))
