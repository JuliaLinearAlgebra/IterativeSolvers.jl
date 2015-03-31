using IterativeSolvers
using Base.Test

#A rank 2 matrix
let A = reshape(1:64,8,8)
    S = rsvd_fnkz(A, 2, ϵ=1e-9)
    @test vecnorm(A - S[:U]*Diagonal(S[:S])*S[:Vt]) < 1e-9
end

#A linearly dependent matrix
let A = [[1:8] [1:8] [1:8] [1:8]]
    S = rsvd_fnkz(A, 4)
    @assert vecnorm(A - S[:U]*Diagonal(S[:S])*S[:Vt]) < 1e-9
end

#SVD of a matrix of zeros should be empty
let
    S = rsvd_fnkz(zeros(10,10), 1)
    @assert S[:U] == zeros(10, 0)
    @assert S[:S] == []
    @assert S[:Vt] == zeros(0, 10)
end
