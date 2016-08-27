using IterativeSolvers
using FactCheck
using Plots

unicodeplots() 	#Use UnicodePlots backend because it is the lightest

facts("ConvergenceHistory") do
    #Just check it doesn't crash
    context("Recipe") do
        A = lu(rand(10, 10))[1]
        b = rand(10)

        for solver in [cg, gmres, lsqr, lsmr, idrs, jacobi, gauss_seidel]
            plot(solver(A,b; log=true)[2])
            @fact true --> true
        end

        for solver in [sor, ssor]
            plot(solver(A,b,1; log=true)[2])
            @fact true --> true
        end

        plot(chebyshev(A,b,1,2; log=true)[2])
        @fact true --> true

        plot(powm(A; log=true)[2])
        @fact true --> true

        plot(invpowm(A; log=true)[2])
        @fact true --> true

		plot(rqi(A; log=true)[2])
        @fact true --> true

        plot(svdl(A; log=true)[3])
        @fact true --> true
    end
end
