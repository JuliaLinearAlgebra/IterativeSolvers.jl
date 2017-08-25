using IterativeSolvers
using Base.Test

@testset "SVD Lanczos" begin

srand(1234567)

#Thick restart methods
@testset "Thick restart with method=$method" for method in (:ritz, :harmonic)
    for T in (Float32, Float64)
        @testset "Diagonal Matrix{$T}" begin
            n = 30
            ns = 5
            tol = 1e-5

            A = full(Diagonal(T[1.0 : n;]))
            q = ones(T, n) / √n
            σ, L, history = svdl(A, nsv=ns, v0=q, tol=tol, reltol=tol, maxiter=n, method=method, vecs=:none, log=true)
            @test isa(history, ConvergenceHistory)
            @test norm(σ - [n : -1.0 : n - 4;]) < 5^2 * 1e-5
            @test_throws ArgumentError svdl(A, nsv=ns, v0=q, tol=tol, reltol=tol, maxiter=n, method=:fakemethod, vecs=:none)

            #Check the singular vectors also
            Σ, L = svdl(A, nsv=ns, v0=q, tol=tol, reltol=tol, maxiter=n, method=method, vecs=:both)

            #The vectors should have the structure
            # [ 0  0 ...  0 ]
            # [    ...      ]
            # [ 0  0 ... ±1 ]
            # [ 0   ...     ]
            # [ 0 ±1 ...  0 ]
            # [±1  0 ...  0 ]
            # and furthermore the signs should be aligned across Σ[:U] and Σ[:V]
            for i = 1 : 5
                Σ[:U][end + 1 - i, i] -= sign(Σ[:U][end + 1 - i, i])
            end
            @test vecnorm(Σ[:U]) < σ[1] * √tol
            for i = 1 : 5
                Σ[:Vt][i, end + 1 - i] -= sign(Σ[:Vt][i, end + 1 - i])
            end
            @test vecnorm(Σ[:U]) < σ[1] * √tol
            @test norm(σ - Σ[:S]) < 2max(tol * ns * σ[1], tol)

            #Issue #55
            let
                σ1, _ = svdl(A, nsv=1, tol=tol, reltol=tol)
                @test abs(σ[1] - σ1[1]) < 2max(tol * σ[1], tol)
            end
        end

        @testset "Rectangular Matrix{$T}" begin
            srand(1)
            m = 300
            n = 200
            k = 5
            l = 10

            A = randn(T, m, n)
            q = randn(T, n) |> x -> x / norm(x)
            σ, L = svdl(A, nsv=k, k=l, v0=q, tol=1e-5, maxiter=30, method=method)
            @test norm(σ - svdvals(A)[1 : k]) < k^2 * 1e-5
        end
    end
end
end #svdl

@testset "BrokenArrowBidiagonal" begin
    B = IterativeSolvers.BrokenArrowBidiagonal([1, 2, 3], [1, 2], Int[])
    @test full(B) == [1 0 1; 0 2 2; 0 0 3]
    @test B[3,3] == 3
    @test B[2,3] == 2
    @test B[3,2] == 0
    @test B[1,3] == 1
    @test size(B) == (3,3)
    @test_throws ArgumentError size(B,3)
    @test_throws BoundsError B[1,5]
end
