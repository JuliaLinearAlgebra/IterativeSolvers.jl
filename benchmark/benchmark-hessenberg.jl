module HessenbergBench

using IterativeSolvers
using BenchmarkTools

function backslash_versus_givens(; n = 1_000, ms = 10 : 10 : 100)
    results = BenchmarkGroup()
    results["backslash"] = BenchmarkGroup()
    results["givens_qr"] = BenchmarkGroup()

    # Example matrix
    A = Tridiagonal(fill(-1.5, n - 1), fill(3.0, n), fill(-1.0, n - 1))

    # Different sizes of the Hessenberg matrix
    for m = ms
        println(m)

        H = zeros(m + 1, m)
        
        # Create an orthonormal basis for the Krylov subspace
        V = rand(n, m + 1)
        V[:, 1] /= norm(V[:, 1])
        
        for i = 1 : m
            # Expand
            V[:, i + 1] = A * V[:, i]
            
            # Orthogonalize
            H[1 : i, i] = view(V, :, 1 : i)' * V[:, i + 1]
            V[:, i + 1] -= view(V, :, 1 : i) * H[1 : i, i]
            
            # Re-orthogonalize
            update = view(V, :, 1 : i)' * V[:, i + 1]
            H[1 : i, i] += update
            V[:, i + 1] -= view(V, :, 1 : i) * update
            
            # Normalize
            H[i + 1, i] = norm(V[:, i + 1])
            V[:, i + 1] /= H[i + 1, i]
        end

        # Fist standard basis vector as rhs
        rhs = [i == 1 ? 1.0 : 0.0 for i = 1 : size(H, 1)]

        # Run the benchmark
        results["givens_qr"][m] = @benchmark A_ldiv_B!(IterativeSolvers.Hessenberg(myH), myRhs) setup = (myH = copy($H); myRhs = copy($rhs))
        results["backslash"][m] = @benchmark $H \ $rhs
    end

    return results
end

end