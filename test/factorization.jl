using IterativeSolvers
using FactCheck

facts("idfact") do
    srand(1)
    M = randn(4,5)
    k = 3

    F = idfact(M, k, 3)
    @fact vecnorm(F.B*F.P-M) --> less_than_or_equal(2svdvals(M)[k+1])
end
