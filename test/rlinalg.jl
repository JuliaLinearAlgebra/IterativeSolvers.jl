using Base.Test
using IterativeSolvers

srand(1)
m = n = 100 
B = randn(m, n)
nB = norm(B)
p = 1e-10 #Probability of failure

for j=1:n
    @test rnorm(B, 2j, p) >= nB
    @test rnorms(B, j, p) >= nB
end

