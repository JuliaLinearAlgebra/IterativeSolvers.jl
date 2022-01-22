VERSION < v"0.7.0-beta2.199" && __precompile__()

"""
Main module for `IterativeSolvers.jl` -- a Julia package for solving linear systems,
eigensystems, and singular value problems.
"""
module IterativeSolvers

include("common.jl")
include("orthogonalize.jl")
include("history.jl")

# Factorizations
include("hessenberg.jl")

# Krylov subspace methods
include("lal.jl")

# Linear solvers
include("stationary.jl")
include("stationary_sparse.jl")
include("cg.jl")
include("minres.jl")
include("bicgstabl.jl")
include("gmres.jl")
include("chebyshev.jl")
include("idrs.jl")
include("qmr.jl")

# Eigensolvers
include("simple.jl")
include("lobpcg.jl")

# SVD solvers
include("svdl.jl")

# Least-squares
include("lsqr.jl")
include("lsmr.jl")

end
