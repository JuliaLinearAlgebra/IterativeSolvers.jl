using Base.Test
using IterativeSolvers
srand(1)
M = randn(4,5)
k = 3

F = idfact(M, k, 3)
vecnorm(F.B*F.P-M) <= 2svdvals(M)[k+1]

