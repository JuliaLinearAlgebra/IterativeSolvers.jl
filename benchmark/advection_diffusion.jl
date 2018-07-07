f(x, y, z) = exp.(x .* y .* z) .* sin.(π .* x) .* sin.(π .* y) .* sin.(π .* z)

function advection_dominated(;N = 50, β = 1000.0)
    # Problem: Δu + βuₓ = f
    # u = 0 on the boundaries
    # f(x, y, z) = exp(xyz) sin(πx) sin(πy) sin(πz)
    # 2nd order central differences (shows serious wiggles)

    # Total number of unknowns
    n = N^3

    # Mesh width
    h = 1.0 / (N + 1)

    # Interior points only
    xs = linspace(0, 1, N + 2)[2 : N + 1]

    # The Laplacian
    Δ = laplace_matrix(Float64, N, 3) ./ -h^2

    # And the dx bit.
    ∂x_1d = spdiagm(-1 => fill(-β / 2h, N - 1), 1 => fill(β / 2h, N - 1))
    ∂x = kron(speye(N^2), ∂x_1d)

    # Final matrix and rhs.
    A = Δ + ∂x
    b = reshape([f(x, y, z) for x ∈ xs, y ∈ xs, z ∈ xs], n)

    A, b
end

function laplace_matrix(::Type{T}, n, dims) where T
    D = second_order_central_diff(T, n)
    A = copy(D)

    for idx = 2 : dims
        A = kron(A, speye(n)) + kron(speye(size(A, 1)), D)
    end

    A
end

second_order_central_diff(::Type{T}, dim) where {T} = convert(
    SparseMatrixCSC{T, Int},
    SymTridiagonal(fill(2 * one(T), dim), fill(-one(T), dim - 1))
)
