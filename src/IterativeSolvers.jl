module IterativeSolvers
include("errors.jl")
include("krylov.jl")

#Linear solvers
include("stationary.jl")
include("cg.jl")
include("gmres.jl")

#Eigensolvers
include("simple.jl")
include("lanczos.jl")

end

