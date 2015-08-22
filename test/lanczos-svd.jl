using IterativeSolvers

#Simple test
let A = full(Diagonal([10.0, 9, 8, 6, 1]))
    @assert norm(svdvals_gkl(A)[1] - svdvals(A)) ≤ 1e-10
end

#Find top singular values of some random triangular matrix
let
    n = 500
    σth = √eps()
    nvals = 6

    A = UpperTriangular(rand(n, n))

    svals = svdvals_gkl(A, nvals)[1]
    svals2 = svdvals(A)[1:nvals]

    @assert norm(svals - svals2) ≤ nvals*σth
end

