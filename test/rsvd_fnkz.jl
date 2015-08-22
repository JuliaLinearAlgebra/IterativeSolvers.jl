using IterativeSolvers
using FactCheck

facts("rsvd_fnkz") do

context("Rank 2 matrix") do
    A = reshape(1:64,8,8)
    S = rsvd_fnkz(A, 2, ϵ=1e-9)
    @fact vecnorm(A - S[:U]*Diagonal(S[:S])*S[:Vt]) --> less_than(1e-9)
end

context("Linearly dependent matrix") do
    A = [collect(1:8) collect(1:8) collect(1:8) collect(1:8)]
    S = rsvd_fnkz(A, 4)
    @fact vecnorm(A - S[:U]*Diagonal(S[:S])*S[:Vt]) --> less_than(1e-9)
end

#SVD of a matrix of zeros should be empty
context("Zero matrix") do
    S = rsvd_fnkz(zeros(10,10), 1)
    @fact S[:U] --> zeros(10, 0)
    @fact S[:S] --> []
    @fact S[:Vt] --> zeros(0, 10)
end

end

