using IterativeSolvers

#Common functions and data structures
include("common.jl")
include("orthogonalize.jl")

# Hessenberg problem solver
include("hessenberg.jl")

#Stationary solvers
include("stationary.jl")

#Conjugate gradients
include("cg.jl")

#BiCGStab(l)
include("bicgstabl.jl")

#MINRES
include("minres.jl")

#GMRES
include("gmres.jl")

#IDRS
include("idrs.jl")

# QMR
include("qmr.jl")

#Chebyshev
include("chebyshev.jl")

#Simple Eigensolvers
include("simple_eigensolvers.jl")
include("lobpcg.jl")

#Golub-Kahan-Lanczos singular values computation
include("svdl.jl")

include("lsqr.jl")
include("lsmr.jl")

#History data structure
include("history.jl")

#Expensive tests - don't run by default
#include("matrixmarket.jl")
#include("matrixcollection.jl")
