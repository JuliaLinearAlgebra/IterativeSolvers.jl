using IterativeSolvers; const IS = IterativeSolvers
using LinearAlgebra
using SparseArrays
using Test

# Equation references and identities from:
# Freund, R. W., & Nachtigal, N. M. (1994). An Implementation of the QMR Method Based on Coupled Two-Term Recurrences. SIAM Journal on Scientific Computing, 15(2), 313–337. https://doi.org/10.1137/0915022

function test_intermediate_identities(ld)
    # 2.7, 3.23
    @test transpose(ld.W[:, 1:end-1]) * ld.V[:, 1:end-1] ≈ ld.D
    # 3.14, 3.26
    @test transpose(ld.Q) * ld.A * ld.P ≈ ld.E
    # 3.10
    @test ld.A * ld.V[:, 1:end-1] ≈ ld.V * ld.L * ld.U
    # 3.7
    @test ld.V[:, 1:end-1] ≈ ld.P * ld.U
    # 3.7
    @test ld.A * ld.P ≈ ld.V * ld.L
    # 5.1
    @test ld.G ≈ ld.U * ld.L[1:end-1, 1:end-1]
    # 3.11
    @test ld.H ≈ ld.L * ld.U

    Γ = Diagonal(ld.γ[1:end-1])
    # 3.15
    @test ld.D * Γ ≈ transpose(ld.D * Γ)
    # 3.16
    @test ld.E * Γ ≈ transpose(ld.E * Γ)

    # 3.8
    @test ld.W[:, 1:end-1] ≈ ld.Q * ((Γ \ ld.U) * Γ)
    # 3.8
    @test ld.At * ld.Q[:, 1:end-1] ≈ ld.W[:, 1:end-1] * (Γ \ ld.L[1:end-1, 1:end-1]) * Γ[1:end-1, 1:end-1]
    # 4.1
    @test ld.A * ld.P[:, 1:end-1] ≈ ld.P * ld.U * ld.L[1:end-1, 1:end-1]
    @test ld.At * ld.Q[:, 1:end-1] ≈ ld.Q * (Γ \ ld.U) * ld.L[1:end-1, 1:end-1] * Γ[1:end-1, 1:end-1]

    # 3.9
    @test all(diag(ld.U) .≈ 1)

    # Definition, by 5.1
    F = transpose(ld.W[:, 1:end-1]) * ld.A * ld.P
    F̃ = transpose(ld.Q) * (A * ld.V[:, 1:end-1])
    # Lemma 5.1
    @test F * Γ ≈ transpose(F̃ * Γ)
    @test ld.F[:, 1:end-1] ≈ F[:, 1:end-1]

    @test all([norm(ld.V[:, i]) for i in axes(ld.V, 2)] .≈ 1)
    @test all([norm(ld.W[:, i]) for i in axes(ld.W, 2)] .≈ 1)
end

function test_finished_identities(ld)
    # 2.7, 3.23
    @test transpose(ld.W) * ld.V ≈ ld.D
    # 3.14, 3.26
    @test transpose(ld.Q[:, 1:end-1]) * ld.A * ld.P[:, 1:end-1] ≈ ld.E
    # 3.10
    @test ld.A * ld.V[:, 1:end-1] ≈ ld.V * ld.L * ld.U[1:end-1, 1:end-1]
    # 3.7
    @test ld.V ≈ ld.P * ld.U
    # 3.7
    @test ld.A * ld.P[:, 1:end-1] ≈ ld.V * ld.L
    # 5.1
    @test ld.G ≈ ld.U * ld.L[1:end-1, 1:end-1]
    # 3.11
    @test ld.H ≈ ld.L * ld.U

    Γ = Diagonal(ld.γ)
    # 3.15
    @test ld.D * Γ ≈ transpose(ld.D * Γ)
    # 3.16
    @test ld.E * Γ[1:end-1, 1:end-1] ≈ transpose(ld.E * Γ[1:end-1, 1:end-1])

    # 3.8
    @test ld.W ≈ ld.Q * ((Γ \ ld.U) * Γ)
    # 3.8
    @test ld.At * ld.Q[:, 1:end-1] ≈ ld.W * (Γ \ ld.L) * Γ[1:end-1, 1:end-1]
    # 4.1
    @test ld.A * ld.P[:, 1:end-1] ≈ ld.P * ld.U * ld.L
    @test ld.At * ld.Q[:, 1:end-1] ≈ ld.Q * (Γ \ ld.U) * ld.L * Γ[1:end-1, 1:end-1]

    # 3.9
    @test all(diag(ld.U) .≈ 1)

    # Definition, by 5.1
    F = transpose(ld.W) * ld.A * ld.P
    F̃ = transpose(ld.Q) * (A * ld.V)
    # Lemma 5.1
    @test F * Γ ≈ transpose(F̃ * Γ)
    @test ld.F[:, 1:end-1] ≈ F[:, 1:end-1]
    @test ld.F̃lastcol ≈ F̃[1:end-1, end]

    # 3.13
    @test ld.L ≈ (ld.D \ Γ) * transpose(ld.U) * (Γ \ transpose(ld.Q)) * A * ld.P[:, end]

    # 3.35
    @test ld.U[1:end-1, 1:end-1] ≈ (ld.E \ transpose(ld.Q[1:end-1])) * (ld.A * ld.V[:, 1:end-1])
    # 3.42
    # Likewise, similar to 3.13 above
    @test ld.L ≈ (ld.D \ transpose(ld.W)) * (ld.A * ld.P[:, 1:end-1])

    @test all([norm(ld.V[:, i]) for i in axes(ld.V, 2)] .≈ 1)
    @test all([norm(ld.W[:, i]) for i in axes(ld.W, 2)] .≈ 1)
end

function test_regular_identities(ld)
    # If all regular vectors, expect:
    # 3.22
    @test abs(sum(triu(ld.U, 2))) < 1e-8
    # 3.20
    @test abs(sum(triu(ld.L, 1))) < 1e-8
    @test abs(sum(triu(ld.D, 1)) + sum(tril(ld.D, -1))) < 1e-8
    @test abs(sum(triu(ld.E, 1)) + sum(tril(ld.E, -1))) < 1e-8
    @test ld.log.PQ_block_count[1] == ld.n
    @test ld.log.VW_block_count[1] == ld.n - 1
end

@testset "A = I" begin
    # A = I terminates immediately (because p1 = v1 -> v2 = Ap1 - v1 = 0)
    A = Diagonal(fill(1.0, 5))
    v = rand(5)
    v ./= norm(v)
    w = rand(5)
    w ./= norm(w)
    ld = IS.LookAheadLanczosDecomp(A, v, w)
    for l in ld end

    @test ld.n == 1
end