using IterativeSolvers
using FactCheck

facts("Randomized linear algebra") do
    srand(1)
    m = n = 100
    B = randn(m, n)
    nB = norm(B)
    p = 1e-5 #Probability of failure

    for j=1:n
        @fact rnorm(B, 2j, p) --> greater_than_or_equal(nB)
        @fact rnorms(B, j, p) --> greater_than_or_equal(nB)
    end

    A = B*B'
    k = 1
    l, u = reigmin(A, k, p)
    @fact l <= eigmin(A) <= u --> true

    l, u = reigmax(A, k, p)
    @fact l <= eigmax(A) <= u --> true

    l, u = rcond(A, k, p)
    @fact l <= cond(A) <= u --> true

    @fact_throws ArgumentError IterativeSolvers.randnn(Char, m)
    @fact_throws ArgumentError IterativeSolvers.randnn(Char, m, n)

end
