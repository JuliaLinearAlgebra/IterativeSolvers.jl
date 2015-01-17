module IterativeSolvers
include("common.jl")
include("krylov.jl")

#Specialized Factorizations
include("factorization.jl")

#Linear solvers
include("stationary.jl")
include("cg.jl")
include("gmres.jl")
include("chebyshev.jl")

#Eigensolvers
include("simple.jl")
include("lanczos.jl")

#Least-squares
include("lsqr.jl")

#Randomized algorithms
include("rlinalg.jl")
include("rsvd.jl")

end
