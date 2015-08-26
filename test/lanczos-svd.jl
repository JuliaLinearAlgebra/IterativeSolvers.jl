using IterativeSolvers
using FactCheck

facts("svdvals_gkl") do

context("Small diagonal matrix") do
    A = full(Diagonal([10.0, 9, 8, 6, 1]))
    @fact norm(svdvals_gkl(A)[1] - svdvals(A)) --> less_than(1e-10)
end

for T in (Float32, Float64)
    context("Medium random UpperTriangular{$T}") do
        srand(1)
        n = 150
        nvals = 6

        Base.Ac_mul_B!(p::Vector, A::UpperTriangular, u::Vector) = Ac_mul_B!(A, copy!(p, u))
        Base.At_mul_B!(p::Vector, A::UpperTriangular, u::Vector) = At_mul_B!(A, copy!(p, u))

        A = UpperTriangular(convert(Matrix{T}, rand(n, n)))
        σth = 0.1*√eps(T)

        svals, = svdvals_gkl(A, nvals, σth=σth)
        @fact length(svals) --> greater_than_or_equal(nvals)

        svals = svals[1:nvals]
        svals2 = svdvals(A)[1:nvals]
        @fact norm(svals - svals2) --> less_than(nvals*σth) "Disagreement in top $nvals singular values:\nsvdvals_gkl\tsvdvals\n"*repr([svals svals2])
    end
end

for T in (Float32, Float64)
    context("Medium wide rectangular Matrix{$T}") do
        srand(1)
        n = 500
        m = 200
        nvals = 10

        A = convert(Matrix{T}, rand(n, m))
        σth = 0.1*√eps(T)

        svals, = svdvals_gkl(A, nvals, σth=σth)
        @fact length(svals) --> greater_than_or_equal(nvals)

        nvals = 4
        svals = svals[1:nvals]
        svals2 = svdvals(A)[1:nvals]
        @fact norm(svals - svals2) --> less_than(nvals*σth) "Disagreement in top $nvals singular values:\nsvdvals_gkl\tsvdvals\n"*repr([svals svals2])
    end

    context("Medium tall rectangular Matrix{$T}") do
        srand(1)
        n = 200
        m = 500
        nvals = 10

        A = convert(Matrix{T}, rand(n, m))
        σth = 0.1*√eps(T)

        svals, = svdvals_gkl(A, nvals, σth=σth)
        @fact length(svals) --> greater_than_or_equal(nvals)

        nvals = 4
        svals = svals[1:nvals]
        svals2 = svdvals(A)[1:nvals]
        @fact norm(svals - svals2) --> less_than(nvals*σth) "Disagreement in top $nvals singular values:\nsvdvals_gkl\tsvdvals\n"*repr([svals svals2])
    end
e
end


end
