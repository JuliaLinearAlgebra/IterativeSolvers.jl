using IterativeSolvers
using FactCheck

facts("rsvd") do

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
    S2 = rsvd(A, 4, 0)

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

context("rrange_adaptive") do
    A = [1. 2 3; 4 5 6; 7 8 9]
    @fact size(IterativeSolvers.rrange_adaptive(A, 3, 1e-3)) --> (3,2)
end

end
