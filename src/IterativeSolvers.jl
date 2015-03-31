module IterativeSolvers

#Documentation support
#see https://github.com/MichaelHatherly/Docile.jl/issues/64
if VERSION < v"0.4-"
    using Docile
    macro doc_mstr(text)
        Base.triplequoted(text)
    end
    macro doc_str(text)
        text
    end
end

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

#Least-squares
include("lsqr.jl")

#Randomized algorithms
include("rlinalg.jl")
include("rsvd.jl")
include("rsvd_fnkz.jl")

end
