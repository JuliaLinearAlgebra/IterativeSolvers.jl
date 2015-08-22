using IterativeSolvers
using FactCheck

facts("svdvals_gkl") do

context("Small diagonal matrix") do
    A = full(Diagonal([10.0, 9, 8, 6, 1]))
    @fact norm(svdvals_gkl(A)[1] - svdvals(A)) --> less_than(1e-10)
end

context("Medium random square matrix") do
    n = 500
    σth = √eps()
    nvals = 6

    A = UpperTriangular(rand(n, n))

    svals = svdvals_gkl(A, nvals)[1]
    svals2 = svdvals(A)[1:nvals]

    @fact norm(svals - svals2) --> less_than(nvals*σth)
end

end
