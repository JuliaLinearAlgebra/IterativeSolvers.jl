using IterativeSolvers
using FactCheck

facts("svdvals_tr") do
#Thick restart methods
for method in (:thick, :harmonic) context("Thick restart with method=$method") do
  for elty in (Float32, Float64)
    context("Diagonal Matrix{$elty}") do
        n = 30
        A = full(Diagonal(elty[1.0:n;]))
        q = convert(Vector{elty}, ones(n)/√n)
        σ, L = svdvals_tr(A, q, 5, 10, tol=1e-5, maxiter=30, method=method)
        @fact norm(σ - [n:-1.0:n-4;]) --> less_than(5^2*1e-5)
    end

    context("Rectangular Matrix{$elty}") do
        srand(1)
        m = 300
        n = 200
        k = 5
        l = 10

        A = convert(Matrix{elty}, randn(m,n))
        q = convert(Vector{elty}, randn(n))|>x->x/norm(x)
        σ, L = svdvals_tr(A, q, k, l, tol=1e-5, maxiter=30, method=method)
        @fact norm(σ - svdvals(A)[1:k]) --> less_than(k^2*1e-5)
    end
  end
end end
end #svdvals_tr

facts("BrokenArrowBidiagonal") do
    B = IterativeSolvers.BrokenArrowBidiagonal([1, 2, 3], [1, 2], Int[])
    @fact full(B) --> [1 0 1; 0 2 2; 0 0 3]
end

