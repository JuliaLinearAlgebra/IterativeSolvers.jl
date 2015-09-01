module IterativeSolvers

isdefined(:__precompile__) && __precompile__(true)

include("common.jl")
include("krylov.jl")

#Specialized factorizations
include("factorization.jl")

#Linear solvers
include("stationary.jl")
include("cg.jl")
include("gmres.jl")
include("chebyshev.jl")

#Eigensolvers
include("simple.jl")
include("lanczos.jl")

#SVD solvers
include("lanczos-svd.jl")
include("lanczos-svd-tr.jl")

#Least-squares
include("lsqr.jl")

#Randomized algorithms
include("rlinalg.jl")
include("rsvd.jl")
include("rsvd_fnkz.jl")

end
