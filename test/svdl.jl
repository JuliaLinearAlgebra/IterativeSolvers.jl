using IterativeSolvers
using FactCheck

facts("svdl") do
#Thick restart methods
for method in (:ritz, :harmonic) context("Thick restart with method=$method") do
  for elty in (Float32, Float64)
    context("Diagonal Matrix{$elty}") do
        n = 30
        ns= 5
        tol=1e-5

        A = full(Diagonal(elty[1.0:n;]))
        q = convert(Vector{elty}, ones(n)/√n)
        σ, L = svdl(A, ns, v0=q, tol=tol, reltol=tol, maxiter=n, method=method, vecs=:none)
        @fact norm(σ - [n:-1.0:n-4;]) --> less_than(5^2*1e-5)
        @fact_throws ArgumentError svdl(A, ns, v0=q, tol=tol, reltol=tol, maxiter=n, method=:fakemethod, vecs=:none)

        #Check the singular vectors also
        Σ, L = svdl(A, ns, v0=q, tol=tol, reltol=tol, maxiter=n, method=method, vecs=:both)

        #The vectors should have the structure
        # [ 0  0 ...  0 ]
        # [    ...      ]
        # [ 0  0 ... ±1 ]
        # [ 0   ...     ]
        # [ 0 ±1 ...  0 ]
        # [±1  0 ...  0 ]
        # and furthermore the signs should be aligned across Σ[:U] and Σ[:V]
        signs = elty[]
        for i=1:5
            Σ[:U][end+1-i, i] -= sign(Σ[:U][end+1-i, i])
        end
        @fact vecnorm(Σ[:U]) --> less_than(σ[1]*√tol)
        for i=1:5
            Σ[:Vt][i, end+1-i] -= sign(Σ[:Vt][i, end+1-i])
        end
        @fact vecnorm(Σ[:U]) --> less_than(σ[1]*√tol)
        @fact norm(σ - Σ[:S]) --> less_than(2max(tol*ns*σ[1], tol))

        #Issue #55
        let
            σ1, _ = svdl(A, 1, tol=tol, reltol=tol)
            @fact abs(σ[1] - σ1[1]) --> less_than(2max(tol*σ[1], tol))
        end
    end

    context("Rectangular Matrix{$elty}") do
        srand(1)
        m = 300
        n = 200
        k = 5
        l = 10

        A = convert(Matrix{elty}, randn(m,n))
        q = convert(Vector{elty}, randn(n))|>x->x/norm(x)
        σ, L = svdl(A, k, k=l, v0=q, tol=1e-5, maxiter=30, method=method)
        @fact norm(σ - svdvals(A)[1:k]) --> less_than(k^2*1e-5)
    end
  end
end end
end #svdl

facts("BrokenArrowBidiagonal") do
    B = IterativeSolvers.BrokenArrowBidiagonal([1, 2, 3], [1, 2], Int[])
    @fact full(B) --> [1 0 1; 0 2 2; 0 0 3]
    @fact B[3,3] --> 3
    @fact B[2,3] --> 2
    @fact B[3,2] --> 0
    @fact B[1,3] --> 1 
    @fact size(B) --> (3,3)
    @fact_throws ArgumentError size(B,3)
    @fact_throws BoundsError B[1,5]
end
