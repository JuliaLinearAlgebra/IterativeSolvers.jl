using IterativeSolvers
using FactCheck
using LinearMaps

srand(1234321)

include("advection_diffusion.jl")

facts("bicgstab(l)") do

for T in (Float32, Float64, Complex64, Complex128)
    context("Matrix{$T}") do

        n = 20
        A = rand(T, n, n) + 15 * eye(T, n)
        x = ones(T, n)
        b = A * x

        for l = (2, 4)
            context("BiCGStab($l)") do

            # Solve without preconditioner
            x1, his1 = bicgstabl(A, b, l, max_mv_products = 100, log = true)
            @fact norm(A * x1 - b) / norm(b) --> less_than(√eps(real(one(T))))

            # Do an exact LU decomp of a nearby matrix
            F = lufact(A + rand(T, n, n))
            x2, his2 = bicgstabl(A, b, Pl = F, l, max_mv_products = 100, log = true)
            @fact norm(A * x2 - b) / norm(b) --> less_than(√eps(real(one(T))))
            end
        end
    end
end
end
