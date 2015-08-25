using IterativeSolvers
using FactCheck

facts("svdvals_gkl") do

context("Small diagonal matrix") do
    A = full(Diagonal([10.0, 9, 8, 6, 1]))
    @fact norm(svdvals_gkl(A)[1] - svdvals(A)) --> less_than(1e-10)
end

context("Medium random square matrix") do
    srand(1)
    n = 500
    σth = 0.1*√eps()
    nvals = 10

    Base.Ac_mul_B!(p::Vector, A::UpperTriangular, u::Vector) = Ac_mul_B!(A, copy!(p, u))
    Base.At_mul_B!(p::Vector, A::UpperTriangular, u::Vector) = At_mul_B!(A, copy!(p, u))

    A = UpperTriangular(rand(n, n))

    svals, = svdvals_gkl(A, nvals, σth=σth)
    @fact length(svals) --> greater_than_or_equal(nvals)

    nvals = 6
    svals = svals[1:nvals]
    svals2 = svdvals(A)[1:nvals]
    @fact norm(svals - svals2) --> less_than(nvals*σth) "Disagreement in top $nvals singular values:\nsvdvals_gkl\tsvdvals\n"*repr([svals svals2])
end

end
