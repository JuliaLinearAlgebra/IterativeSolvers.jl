
export rsvd_fnkz

immutable OuterProduct{T}
    X :: Matrix{T}
    Y :: Matrix{T}
end

X ⊗ Y = OuterProduct{promote_type(eltype(X), eltype(Y))}(X, Y)

"""
    rsvd_fnkz(A, k)

Compute the randomized SVD by iterative refinement from randomly selected
columns/rows.

# Arguments

`A`: matrix whose SVD is desired.

`k::Int`: desired rank of approximation (`k ≤ min(m, n)`).

## Keywords

`l::Int = k`: number of columns/rows to sample at each iteration (`1 ≤ l ≤ k`).

`N::Int = minimum(size(A))`: maximum number of iterations.

`ϵ::Real = prod(size(A))*eps()`: relative threshold for convergence, as
measured by growth of the spectral norm.

`method::Symbol = :eig`: problem to solve.
* `:eig`: eigenproblem.
* `:svd`: singular problem.

`verbose::Bool = false`: print convergence information at each iteration.

# Output

SVD object of `rank ≤ k`.

# References

```bibtex
@inproceedings{,
    author={Friedland, S. and Niknejad, A. and Kaveh, Mostafa and Zare, H.},
    booktitle={System of Systems Engineering, 2006 IEEE/SMC International Conference on},
    title={Fast Monte-Carlo low rank approximations for matrices},
    year={2006},
    month={April},
    pages={218--223},
    doi={10.1109/SYSOSE.2006.1652299}
}
```
"""
function rsvd_fnkz(A, k::Int;
    l::Int=k, N::Int=minimum(size(A)), verbose::Bool=false,
    method::Symbol=:eig,
    ϵ::Real=prod(size(A))*eps(real(float(one(eltype(A))))))

    const m, n = size(A)
    const dosvd = method == :svd
    @assert 1 ≤ l ≤ k ≤ min(m, n)

    #Initialize k-rank approximation B₀ according to (III.9)
    const tallandskinny = m ≥ n
    B₀ = if tallandskinny
        #Select k columns of A
        π = randperm(n)[1:k]
        A[:, π]
    else
        warn("Not implemented")
        #Select k rows of A
        π = randperm(m)[1:k]
        A[π, :]
    end
    X, RRR = qr(B₀)
    X = X[:, abs(diag(RRR)) .> ϵ] #Remove linear dependent columns
    B = tallandskinny ? X ⊗ A'X : X ⊗ X'A

    oldnrmB = 0.0
    Λ = 0
    X̃ = 0
    V = 0
    for t=1:N
        π = randperm(k)[1:l]
        #Update B using Theorem 2.4
        X, RRR = qr([B.X A[:, π]]) #We are doing more work than needed here
        X = X[:, abs(diag(RRR)) .> ϵ] #Remove linearly dependent columns
        Y = A'X
        if dosvd
            S = svdfact!(Y)
            Λ = S[:S].^2
            B = (X * S[:V]) ⊗ (S[:U]*Diagonal(S[:S]))
        else
            E = eigfact!(Symmetric(Y'Y))
            π = sortperm(E[:values], rev=true)[1:min(length(E[:values]), k)]
            Λ = E[:values][π]
            O = E[:vectors][:, π] #Eq. 2.6
            X̃ = X * O
            B = X̃ ⊗ A'X̃
        end

        length(Λ)==0 && break

        nrmB = √Λ[1]
        verbose && println("iteration $t: Norm: $nrmB")
        oldnrmB > (1-ϵ)*nrmB && break
        oldnrmB = nrmB
    end
    for i=1:size(B.Y, 2)
        scale!(view(B.Y, :, i), 1/√Λ[i])
    end
    Base.LinAlg.SVD(B.X, √Λ, B.Y')
end
