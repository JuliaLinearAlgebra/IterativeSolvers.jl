module IterativeSolvers
include("common.jl")
include("krylov.jl")

#Orthogonalization routines
include("qr.jl")

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

end
