__precompile__(true)
"""
Main module for `IterativeSolvers.jl` -- a Julia package for solving linear systems,
eigensystems, and singular value problems.
"""
module IterativeSolvers

using SugarBLAS

include("common.jl")
include("orthogonalize.jl")
include("history.jl")

# Factorizations
include("hessenberg.jl")

# Linear solvers
include("stationary.jl")
include("stationary_sparse.jl")
include("cg.jl")
include("minres.jl")
include("bicgstabl.jl")
include("gmres.jl")
include("chebyshev.jl")
include("idrs.jl")

# Eigensolvers
include("simple.jl")

# SVD solvers
include("svdl.jl")

# Least-squares
include("lsqr.jl")
include("lsmr.jl")

end
