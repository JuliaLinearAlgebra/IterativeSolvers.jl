using IterativeSolvers
using FactCheck

srand(1234321)

for rsvd in [IterativeSolvers.rsvdfact]
    facts(rsvd) do

        context("Small wide rectangular") do
            A = [1. 2 3 4; 5 6 7 8]
            S1 = svdfact(A)
            S2 = rsvd(A, 2, 0)

            @fact vecnorm(abs(S1[:U]) - abs(S2[:U])) --> less_than( √(eps()))
            @fact vecnorm(abs(S1[:Vt]) - abs(S2[:Vt])) --> less_than(√(eps()))
            @fact norm(S1[:S] - S2[:S]) --> less_than(√(eps()))
        end

        context("Small tall rectangular") do
            A = [1. 2; 3 4; 5 6; 7 8]
            S1 = svdfact(A)
            S2 = rsvd(A, 2, 0)

            @fact vecnorm(abs(S1[:U]) - abs(S2[:U])) --> less_than(√(eps()))
            @fact vecnorm(abs(S1[:Vt]) - abs(S2[:Vt])) --> less_than(√(eps()))
            @fact norm(S1[:S] - S2[:S]) --> less_than(√(eps()))
        end

        context("Small square") do
            A = [1. 2 3; 4 5 6; 7 8 9]
            S1 = svdfact(A)
            S2 = rsvd(A, 3, 0)

            @fact vecnorm(abs(S1[:U]) - abs(S2[:U]))  --> less_than(√(eps()))
            @fact vecnorm(abs(S1[:Vt]) - abs(S2[:Vt]))  --> less_than(√(eps()))
            @fact norm(S1[:S] - S2[:S])  --> less_than(√(eps()))
        end

        context("Low rank") do #Issue #42
            n = 10
            r = 2
            A = randn(n, r)*randn(r, n)
            S = svdvals(A)
            for nvals = 1:r
                S1= IterativeSolvers.rsvdvals(A, nvals, r-nvals)
                for i = 1:nvals
                    @fact abs(S[i] - S1[i]) --> less_than(n^2*r*eps())
                end
            end
        end

        context("rrange_adaptive") do
            A = [1. 2 3; 4 5 6; 7 8 9]
            @fact size(IterativeSolvers.rrange_adaptive(A, 3, 1e-3)) --> (3,2)
        end

        context("rrange") do
            A = [1. 2 3; 4 5 6; 7 8 9]
            @fact_throws ArgumentError IterativeSolvers.rrange(A, 20)
        end
    end
end
