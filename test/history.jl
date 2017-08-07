using IterativeSolvers
using Plots

#Just check it doesn't crash
@testset "ConvergenceHistory" begin

srand(1234321)

#Use UnicodePlots backend because it is the lightest
unicodeplots()

A = lu(rand(10, 10))[1]
b = rand(10)

for solver in (cg, gmres, minres, lsqr, lsmr, idrs, jacobi, gauss_seidel)
    plot(solver(A, b; log=true)[2])
end

for solver in (sor, ssor)
    plot(solver(A, b, 1.0; log=true)[2])
end

plot(bicgstabl(A, b, 2, log=true)[2])
plot(chebyshev(A, b, 0.5, 1.5; log=true)[2])
plot(powm(A; log=true)[2])
plot(invpowm(A; log=true)[2])
plot(svdl(A; log=true)[3])
end
