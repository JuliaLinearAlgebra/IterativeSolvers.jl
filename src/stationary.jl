#Stationary iterative methods
#Templates, section 2.2
export jacobi, jacobi!, gauss_seidel, gauss_seidel!, sor, sor!, ssor, ssor!

import Base.LinAlg.SingularException

function check_diag(A::AbstractMatrix)
    for i = 1 : size(A, 1)
        if iszero(A[i,i])
            throw(SingularException(i))
        end
    end
end

####################
# API method calls #
####################

jacobi(A::AbstractMatrix, b; kwargs...) = jacobi!(zerox(A, b), A, b; kwargs...)

function jacobi!(x, A::AbstractMatrix, b; maxiter::Int=size(A, 2))
    check_diag(A)
    iterable = DenseJacobiIterable(A, x, similar(x), b, maxiter)
    for _ = iterable end
    x
end

mutable struct DenseJacobiIterable{matT,vecT}
    A::matT
    x::vecT
    next::vecT
    b::vecT
    maxiter::Int
end

start(::DenseJacobiIterable) = 1
done(it::DenseJacobiIterable, iteration::Int) = iteration > it.maxiter
function next(j::DenseJacobiIterable, iteration::Int)
    n = size(j.A, 1)

    copy!(j.next, j.b)
    
    # Computes next = b - (A - D)x
    for col = 1 : n
        @simd for row = 1 : col - 1
            @inbounds j.next[row] -= j.A[row, col] * j.x[col]
        end

        @simd for row = col + 1 : n
            @inbounds j.next[row] -= j.A[row, col] * j.x[col]
        end
    end

    # Computes x = D \ next
    for col = 1 : n
        @inbounds j.x[col] = j.next[col] / j.A[col, col]
    end

    nothing, iteration + 1
end

####################
# API method calls #
####################

gauss_seidel(A::AbstractMatrix, b; kwargs...) = gauss_seidel!(zerox(A, b), A, b; kwargs...)

function gauss_seidel!(x, A::AbstractMatrix, b; maxiter::Int=size(A,2))
    check_diag(A)
    iterable = DenseGaussSeidelIterable(A, x, b, maxiter)
    for _ = iterable end
    x
end

mutable struct DenseGaussSeidelIterable{matT,vecT}
    A::matT
    x::vecT
    b::vecT
    maxiter::Int
end

start(::DenseGaussSeidelIterable) = 1
done(it::DenseGaussSeidelIterable, iteration::Int) = iteration > it.maxiter

function next(s::DenseGaussSeidelIterable, iteration::Int)
    n = size(s.A, 1)
    
    for col = 1 : n
        @simd for row = 1 : col - 1
            @inbounds s.x[row] -= s.A[row, col] * s.x[col]
        end

        s.x[col] = s.b[col]
    end

    for col = 1 : n
        @inbounds s.x[col] /= s.A[col, col]
        @simd for row = col + 1 : n
            @inbounds s.x[row] -= s.A[row, col] * s.x[col]
        end
    end

    nothing, iteration + 1
end

####################
# API method calls #
####################

sor(A::AbstractMatrix, b, ω::Real; kwargs...) =
    sor!(zerox(A, b), A, b, ω; kwargs...)

function sor!(x, A::AbstractMatrix, b, ω::Real; maxiter=size(A, 2))
    check_diag(A)
    iterable = DenseSORIterable(A, x, similar(x), b, ω, maxiter)
    for _ = iterable end
    x
end

mutable struct DenseSORIterable{matT,vecT,numT}
    A::matT
    x::vecT
    tmp::vecT
    b::vecT
    ω::numT
    maxiter::Int
end

start(::DenseSORIterable) = 1
done(it::DenseSORIterable, iteration::Int) = iteration > it.maxiter
function next(s::DenseSORIterable, iteration::Int)
    n = size(s.A, 1)

    for col = 1 : n
        @simd for row = 1 : col - 1
            @inbounds s.tmp[row] -= s.A[row, col] * s.x[col]
        end

        s.tmp[col] = s.b[col]
    end

    for col = 1 : n
        @inbounds s.x[col] += s.ω * (s.tmp[col] / s.A[col, col] - s.x[col])
        @simd for row = col + 1 : n
            @inbounds s.tmp[row] -= s.A[row, col] * s.x[col]
        end
    end

    nothing, iteration + 1
end

####################
# API method calls #
####################

ssor(A::AbstractMatrix, b, ω::Real; kwargs...) =
    ssor!(zerox(A, b), A, b, ω; kwargs...)

function ssor!(x, A::AbstractMatrix, b, ω::Real; maxiter::Int=size(A,2))
    check_diag(A)
    iterable = DenseSSORIterable(A, x, similar(x), b, ω, maxiter)
    for _ = iterable end
    x
end

mutable struct DenseSSORIterable{matT,vecT,numT}
    A::matT
    x::vecT
    tmp::vecT
    b::vecT
    ω::numT
    maxiter::Int
end

start(::DenseSSORIterable) = 1
done(it::DenseSSORIterable, iteration::Int) = iteration > it.maxiter
function next(s::DenseSSORIterable, iteration::Int)
    n = size(s.A, 1)

    for col = 1 : n
        @simd for row = 1 : col - 1
            @inbounds s.tmp[row] -= s.A[row, col] * s.x[col]
        end

        s.tmp[col] = s.b[col]
    end

    for col = 1 : n
        @inbounds s.x[col] += s.ω * (s.tmp[col] / s.A[col, col] - s.x[col])
        @simd for row = col + 1 : n
            @inbounds s.tmp[row] -= s.A[row, col] * s.x[col]
        end
    end

    for col = n : -1 : 1
        s.tmp[col] = s.b[col]
        @simd for row = col + 1 : n
            @inbounds s.tmp[row] -= s.A[row, col] * s.x[col]
        end
    end

    for col = n : -1 : 1
        @simd for row = 1 : col - 1
            @inbounds s.tmp[row] -= s.A[row, col] * s.x[col]
        end

        @inbounds s.x[col] += s.ω * (s.tmp[col] / s.A[col, col] - s.x[col])
    end

    nothing, iteration + 1
end

#################
# Documentation #
#################

let
#Initialize parameters
doc1_call = """    jacobi(A, b)
"""
doc1!_call = """    jacobi!(x, A, b)
"""
doc2_call = """    gauss_seidel(A, b)
"""
doc2!_call = """    gauss_seidel!(x, A, b)
"""
doc3_call = """    sor(A, b, ω)
"""
doc3!_call = """    sor!(x, A, b, ω)
"""
doc4_call = """    ssor(A, b, ω)
"""
doc4!_call = """    ssor!(x, A, b, ω)
"""
doc1_msg = "Solve A*x=b with the Jacobi method."
doc2_msg = "Solve A*x=b with the Gauss-Seidel method."
doc3_msg = "Solve A*x=b with the successive overrelaxation method."
doc4_msg = "Solve A*x=b with the symmetric successive overrelaxation method."
doc1!_msg = "Overwrite `x`.\n\n" * doc1_msg
doc2!_msg = "Overwrite `x`.\n\n" * doc2_msg
doc3!_msg = "Overwrite `x`.\n\n" * doc3_msg
doc4!_msg = "Overwrite `x`.\n\n" * doc4_msg
doc1_arg = ""
doc2_arg = ""
doc3_arg = "`shift::Number=0`: shift to be applied to matrix A."
doc4_arg = "`shift::Number=0`: shift to be applied to matrix A."
doc1!_arg = "`x`: initial guess, overwrite final estimation."
doc2!_arg = "`x`: initial guess, overwrite final estimation."
doc3!_arg = "`x`: initial guess, overwrite final estimation.\n\n$doc3_arg"
doc4!_arg = "`x`: initial guess, overwrite final estimation.\n\n$doc4_arg"

doc1_version = (jacobi, doc1_call, doc1_msg, doc1_arg)
doc2_version = (gauss_seidel, doc2_call, doc2_msg, doc2_arg)
doc3_version = (sor, doc3_call, doc3_msg, doc3_arg)
doc4_version = (ssor, doc4_call, doc4_msg, doc4_arg)
doc1!_version = (jacobi!, doc1!_call, doc1!_msg, doc1!_arg)
doc2!_version = (gauss_seidel!, doc2!_call, doc2!_msg, doc2!_arg)
doc3!_version = (sor!, doc3!_call, doc3!_msg, doc3!_arg)
doc4!_version = (ssor!, doc4!_call, doc4!_msg, doc4!_arg)

i=0
docstring = Vector(8)

#Build docs
for (func, call, msg, arg) in [doc1_version, doc2_version, doc3_version, doc4_version,
                                doc1!_version, doc2!_version, doc3!_version, doc4!_version]
i+=1
docstring[i] = """
$call

$msg

`ch` is a `ConvergenceHistory` object. Otherwise it will only return `x`.
If `log` is set to `true` is given, method will output a tuple `x, ch`.

# Arguments

$arg

`A`: linear operator.

## Keywords

`tol::Real = size(A,2)^3*eps()`: stopping tolerance.

`maxiter::Integer = size(A,2)^2`: maximum number of iterations.

`verbose::Bool = false`: verbose flag.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

# Output

**if `log` is `false`**

`x`: approximated solution.

**if `log` is `true`**

`x`: approximated solution.

`ch`: convergence history.

**ConvergenceHistory keys**

`:tol` => `::Real`: stopping tolerance.

`:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> jacobi
@doc docstring[5] -> jacobi!
@doc docstring[2] -> gauss_seidel
@doc docstring[6] -> gauss_seidel!
@doc docstring[3] -> sor
@doc docstring[7] -> sor!
@doc docstring[4] -> ssor
@doc docstring[8] -> ssor!
end
