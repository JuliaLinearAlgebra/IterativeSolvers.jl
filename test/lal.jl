using IterativeSolvers; const IS = IterativeSolvers
using LinearAlgebra
using SparseArrays
using Test

# Equation references and identities from:
# Freund, R. W., & Nachtigal, N. M. (1994). An Implementation of the QMR Method Based on Coupled Two-Term Recurrences. SIAM Journal on Scientific Computing, 15(2), 313–337. https://doi.org/10.1137/0915022

function _iterate_and_test_intermediates(ld)
    # iterates through ld and collects the intermediate matrices by appending the last
    # row and column to the matrix being built
    # returns a NamedTuple with the same names as `ld`
    # All intermediates will be equal to the definitions at iteration n
    # Note, at the end of iteration n, the following are set to n+1:
    # ρ, ξ, γ, w, v, w̃, ṽ, w̃tṽ, wtv
    # At the current commit, this seems redundant, but in future commits the
    # calculations will sparsify and only small vectors will be kept
    # Because we don't explicitly construct H and G, we do not attempt collecting them

    P = copy(ld.P) # size (N, n)
    Q = copy(ld.Q) # size (N, n)
    W = copy(ld.W) # size (N, n+1)
    V = copy(ld.V) # size (N, n+1)
    γ = copy(ld.γ) # size n+1
    D = copy(ld.D) # size (n, n)
    E = copy(ld.E) # size (n, n)
    F = copy(ld.F) # size (n, n)
    F̃ = copy(ld.F) # size (n, n)
    U = copy(ld.U) # size (n, n)
    L = copy(ld.L) # size (n+1, n)
    Γ = Diagonal(γ) # size (n+1, n+1)

    for (i, l) in enumerate(ld)
        P = [P ld.p]
        Q = [Q ld.q]
        V = [V ld.v]
        W = [W ld.w]
        push!(γ, ld.γ[end])
        Γ = Diagonal(γ) # size (n+1, n+1)
        if (i==1)
            D = ld.D[:, :]
            E = ld.E[:, :]
            F = ld.F[:, :]
            F̃ = (transpose(F * Γ[1:end-1, 1:end-1]) / Γ[1:end-1, 1:end-1])
            U = ld.U[:, :]
            L = ld.L[:, :]
        else
            D = [D ld.D[1:end-1, end]; ld.D[end:end, :]]
            E = [E ld.E[1:end-1, end]; ld.E[end:end, :]]
            F = [F ld.F[1:end-1, end]; ld.F[end:end, :]]
            # F̃ row is not explicitly calculated, so we calculate from F using Lemma 5.1
            F̃ = [F̃ ld.F̃lastcol; (transpose(F * Γ[1:end-1, 1:end-1]) / Γ[1:end-1, 1:end-1])[end:end, :]]
            U = [U ld.U[1:end-1, end]; ld.U[end:end, :]]
            L = [L ld.L[1:end-1, end]; ld.L[end:end, :]]

            results = (P=P, Q=Q, W=W, V=V, γ=γ, Γ=Γ, D=D, E=E, F=F, F̃=F̃, U=U, L=L, n=ld.n, A=ld.A, At=ld.At)
            _test_lal_identities(results)
        end
    end

    return (P=P, Q=Q, W=W, V=V, γ=γ, Γ=Γ, D=D, E=E, F=F, F̃=F̃, U=U, L=L, n=ld.n, A=ld.A, At=ld.At)
end

function _test_lal_identities(ld)
    # Note: V, W, Γ are all at n+1, but the values at n are use for some identities
    Vn = ld.V[:, 1:end-1]
    Wn = ld.W[:, 1:end-1]
    Γn = ld.Γ[1:end-1, 1:end-1]
    E_singular = cond(ld.E) > 1e12 # ad-hoc singularity definition
    D_singular = cond(ld.D) > 1e12
    # 2.7, 3.23
    @test transpose(Wn) * Vn ≈ ld.D
    # 3.7
    @test Vn ≈ ld.P * ld.U
    # 3.7
    @test ld.A * ld.P ≈ ld.V * ld.L
    # 3.8
    @test Wn ≈ ld.Q * ((Γn \ ld.U) * Γn)
    # 3.8
    @test ld.At * ld.Q ≈ ld.W * (ld.Γ \ ld.L) * Γn
    # 3.9
    @test all(diag(ld.U) .≈ 1)
    # 3.10
    @test ld.A * Vn ≈ ld.V * ld.L * ld.U
    # 3.13
    if !D_singular
        @test ld.L[1:end-1, end] ≈ (ld.D \ Γn) * transpose(ld.U) * (Γn \ transpose(ld.Q)) * ld.A * ld.P[:, end]
    end
    # 3.14, 3.26
    @test ld.E ≈ transpose(ld.Q) * ld.A * ld.P

    # 3.15
    @test ld.D * Γn ≈ transpose(ld.D * Γn)
    # 3.16
    @test ld.E * Γn ≈ transpose(ld.E * Γn)


    # 3.35
    if !E_singular
        @test ld.U[:, end] ≈ (ld.E \ transpose(ld.Q)) * (ld.A * Vn[:, end])
    end
    # 3.42
    # Likewise, similar to 3.13 above
    if !D_singular
        @test ld.L[1:end-1, end] ≈ (ld.D \ transpose(Wn)) * (ld.A * ld.P[:, end])
    end

    # 4.1
    @test ld.A * ld.P[:, 1:end-1] ≈ ld.P * ld.U * ld.L[1:end-1, 1:end-1]
    @test ld.At * ld.Q[:, 1:end-1] ≈ ld.Q * (Γn \ ld.U) * ld.L[1:end-1, 1:end-1] * Γn[1:end-1, 1:end-1]

    # Definition, by 5.1
    F = transpose(Wn) * ld.A * ld.P
    F̃ = transpose(ld.Q) * (ld.A * Vn)
    # Lemma 5.1
    @test F * Γn ≈ transpose(F̃ * Γn)
    @test ld.F ≈ F
    @test ld.F̃ ≈ F̃

    @test all([norm(ld.V[:, i]) for i in axes(ld.V, 2)] .≈ 1)
    @test all([norm(ld.W[:, i]) for i in axes(ld.W, 2)] .≈ 1)
end

function test_regular_identities(ld, log; early_exit=false)
    # If all regular vectors, expect:
    # Eq. 3.22, U is upper triangular (bi-diagonal)
    @test tril(ld.U, 1) ≈ triu(ld.U)
    # Eq. 3.20, L is upper Hessenburg with nothing above diagonal
    @test triu(ld.L, -1) ≈ tril(ld.L)
    @test Diagonal(diag(ld.D)) ≈ ld.D
    @test Diagonal(diag(ld.E)) ≈ ld.E
    @test log.PQ_block_count[1] == ld.n
    if early_exit
        @test log.VW_block_count[1] == ld.n-1
    else
        @test log.VW_block_count[1] == ld.n
    end
end

@testset "A = I" begin
    # A = I terminates immediately (because p1 = v1 -> v2 = Ap1 - v1 = 0)
    A = Diagonal(fill(1.0, 5))
    v = rand(5)
    w = rand(5)
    ld = IS.LookAheadLanczosDecomp(A, v, w; log=true, vw_normalized=false)
    ld_results = _iterate_and_test_intermediates(ld)

    @test ld.n == 1
    @test ld_results.V ≈ v
    @test ld_results.W ≈ w
    @test isempty(ld_results.P)
    @test isempty(ld_results.Q)

    test_regular_identities(ld_results, ld.log; early_exit=true)
end

@testset "Dense Hermitian Matrix" begin
    # Hermitian matrices avoid breakdown in exact arithmetic
    A = [i == j ? 20.0+0.0im : i < j ? -1.0im : 1.0im for i in 1:10, j in 1:10]
    v = fill(0.0im, 10)
    v[1] = 1.0
    w = fill(0.0im, 10)
    w[1:2] .= 1.0
    ld = IS.LookAheadLanczosDecomp(A, v, w; log=true, vw_normalized=false)
    ld_results = _iterate_and_test_intermediates(ld)
    test_regular_identities(ld_results, ld.log)
end
    

@testset "Sparse Tri-diagonal Matrix" begin
    A = Tridiagonal(-ones(9), 2*ones(10), -ones(9))
    v = fill(0.0im, 10)
    v[1] = 1.0
    w = fill(0.0im, 10)
    w[1] = 1.0
    ld = IS.LookAheadLanczosDecomp(A, v, w; log=true, vw_normalized=false)
    ld_results = _iterate_and_test_intermediates(ld)
    test_regular_identities(ld_results, ld.log)
end

@testset "Cyclic Circulant Matrix" begin
    # 4-cyclic circulant matrix, creates blocks
    I1 = Diagonal(fill(1.0, 3))
    I2 = Diagonal(fill(1.0, 7))
    I3 = Diagonal(fill(1.0, 3))
    I4 = Diagonal(fill(1.0, 5))
    B1 = rand(size(I1, 1), size(I4, 2))
    B2 = rand(size(I2, 1), size(I1, 2))
    B3 = rand(size(I3, 1), size(I2, 2))
    B4 = rand(size(I4, 1), size(I3, 2))
    Z12 = fill(0.0, size(I1, 1), size(I2, 2))
    Z13 = fill(0.0, size(I1, 1), size(I3, 2))
    Z14 = fill(0.0, size(I1, 1), size(I4, 2))
    Z23 = fill(0.0, size(I2, 1), size(I3, 2))
    Z24 = fill(0.0, size(I2, 1), size(I4, 2))
    Z34 = fill(0.0, size(I3, 1), size(I4, 2))
    Ac = [
        I1   Z12  Z13 B1
        B2   I2   Z23 Z24
        Z13' B3   I3  Z34
        Z14' Z24' B4  I4
    ]
    v = [ones(size(I1, 1)); fill(0.0, size(Ac, 1) - size(I1, 1))]
    w = [ones(size(I1, 1)); fill(0.0, size(Ac, 1) - size(I1, 1))]
    ld = IS.LookAheadLanczosDecomp(Ac, v, w; log=true, vw_normalized=false, max_block_size=8, verbose=true)
    ld_results = _iterate_and_test_intermediates(ld)
end

