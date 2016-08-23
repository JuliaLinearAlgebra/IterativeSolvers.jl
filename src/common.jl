import  Base: eltype, empty!, eps, length, ndims, push!, real, size, *, \,
        A_mul_B!, Ac_mul_B, Ac_mul_B!, getindex, setindex!, push!

export  A_mul_B

# Improve readability of iterative methods
\(f::Function, b) = f(b)
*(f::Function, b) = f(b)

#### Type-handling
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))
Amultype(A, x) = typeof(one(eltype(A))*one(eltype(x)))
if VERSION < v"0.4.0-dev+6068"
    real{T<:Real}(::Type{Complex{T}}) = T
    real{T<:Real}(::Type{T}) = T
end
eps{T<:Real}(::Type{Complex{T}}) = eps(T)

function randx(A, b)
    T = Adivtype(A, b)
    x = initrand!(Array(T, size(A, 2)))
end

function zerox(A, b)
    T = Adivtype(A, b)
    x = zeros(T, size(A, 2))
end

#### Numerics
function update!(x, α::Number, p::AbstractVector)
    for i = 1:length(x)
        x[i] += α*p[i]
    end
    x
end

function initrand!(v::Vector)
    _randn!(v)
    nv = norm(v)
    for i = 1:length(v)
        v[i] /= nv
    end
    v
end
_randn!(v::Array{Float64}) = randn!(v)
_randn!(v) = copy!(v, randn(length(v)))

#### Reporting
"""
Store general and in-depth information about an iterative method.

**Fields**

* `mvps::Int`: number of matrix vector products.
* `mtvps::Int`: number of transposed matrix-vector products
* `iters::Int`: iterations taken by the method.
* `restart::T`: restart relevant information.
    - `T == Int`: iterations per restart.
    - `T == Void`: methods without restarts.

* `isconverged::Bool`: convergence of the method.
* `tol::Dict{Symbol, Real}`: tolerances of the method.
* `data::Dict{Symbol, VecOrMat}`: iteration information of a method. It usually
contains residuals, but can have other information, e.g. ritz values in [svdl](@ref).

**Constructors**

    ConvergenceHistory()
    ConvergenceHistory(restart)

Create `ConvergenceHistory` with empty fields.

**Arguments**

* `restart`: number of iterations per restart.

**Plots**

Supports plots using the `Plots.jl` package via a type recipe. Vectors are
ploted as series and matrices as scatterplots.

**Implements**

* `Base`: `getindex`, `setindex!`, `push!`

"""
type ConvergenceHistory{T,K}
    mvps::Int
    mtvps::Int
    iters::Int
    restart::K
    isconverged::Bool
    tol::Dict{Symbol, Real}
    data::Dict{Symbol, VecOrMat}
end
function ConvergenceHistory(;restart=nothing, partial=true)
    ConvergenceHistory{!partial,typeof(restart)}(0,0,0,restart,false,
        Dict{Symbol, Real}(), Dict{Symbol, VecOrMat}()
        )
end

"""
Stores information of the current iteration.
"""
typealias PartialHistory ConvergenceHistory{false}

"""
Stores the information of all the iterations.
"""
typealias CompleteHistory ConvergenceHistory{true}

"""
History without resets.
"""
typealias UnrestartedHistory{T} ConvergenceHistory{T, Void}

"""
History with resets.
"""
typealias RestartedHistory{T} ConvergenceHistory{T, Int}

"""
    getindex(ch, s)

Get collection or tolerance associated with key `s` in `ch::ConvergenceHistory`.

    getindex(ch, s, kwargs...)

Access elements of the collection associated with key `s` in `ch::ConvergenceHistory`.
"""
function getindex(ch::PartialHistory, s::Symbol)
    if haskey(ch.tol, s)
        ch.tol[s]
    else
        ch.data[s][1]
    end
end
function getindex(ch::CompleteHistory, s::Symbol)
    if haskey(ch.tol, s)
        ch.tol[s]
    else
        ch.data[s]
    end
end
getindex(ch::CompleteHistory, s::Symbol, kwargs...) = ch.data[s][kwargs...]

"""
    setindex!(ch, tol, s)

Set tolerance value associated with `s` in `ch::ConvergenceHistory` to `tol`.

    setindex!(ch, val, s, kwargs...)

Set collection element associated with key `s` in `ch::ConvergenceHistory` to val.
"""
setindex!(ch::ConvergenceHistory, val::Real, s::Symbol) = ch.tol[s] = val
setindex!(ch::CompleteHistory, val, s::Symbol, kwargs...) = ch.data[s][kwargs...] = val

"""
    push!(ch, key, val)

Push contents of `val` to collection associated with `key` in `ch::ConvergenceHistory`.
"""
function push!(ch::ConvergenceHistory, key::Symbol, vec::Union{Vector,Tuple})
    matrix = ch.data[key]
    width = size(matrix,2)
    iter = isa(ch,CompleteHistory) ? ch.iters : 1
    base = (iter-1)*width
    for i in 1:min(width,length(vec))
        matrix[base+i] = vec[i]
    end
end
push!(ch::ConvergenceHistory, key::Symbol, vec) = push_custom_data!(ch, key, vec)

push_custom_data!(ch::PartialHistory, key::Symbol, val) = ch.data[key][1] = val
push_custom_data!(ch::CompleteHistory, key::Symbol, val) = ch.data[key][ch.iters] = val

"""
    reserve!(ch, key, maxiter)
    reserve!(typ, key, maxiter)
    reserve!(ch, key, maxiter, size)
    reserve!(typ, ch, key, maxiter, size)

Reserve space for per iteration data in `ch`. If size is provided, intead of a
vector it will reserve matrix of dimensions `(maxiter, size)`.

**Arguments**

* `typ::Type`: Type of the elements to store. Defaults to `Float64` when not given.
* `ch::ConvergenceHistory`: convergence history.
* `key::Symbol`: key used to identify the data.
* `maxiter::Int`: number of iterations to save space for.
* `size::Int`: number of elements to store with the `key` identifier.

"""
function reserve!(ch::ConvergenceHistory, key::Symbol, maxiter::Int, kwargs...)
    len = isa(ch,CompleteHistory) ? maxiter : 1
    _reserve!(Float64, ch, key, len, kwargs...)
end
function reserve!(typ::Type, ch::ConvergenceHistory, key::Symbol, maxiter::Int, kwargs...)
    len = isa(ch,CompleteHistory) ? maxiter : 1
    _reserve!(typ, ch, key, len, kwargs...)
end

function _reserve!(typ::Type, ch::ConvergenceHistory, key::Symbol, len::Int)
    ch.data[key] = Vector{typ}(len)
end
function _reserve!(typ::Type, ch::ConvergenceHistory, key::Symbol, len::Int, size::Int)
    ch.data[key] = Matrix{typ}(len,size)
end

"""
    shrink!(ml)

shrinks the reserved space for `MethodLog` `ch` to the space actually used to log.
"""
shrink!(::PartialHistory) = nothing
function shrink!(ch::CompleteHistory)
    for key in datakeys(ch)
        elem = ch.data[key]
        if isa(elem, Vector)
            resize!(elem, ch.iters)
        elseif isa(elem, Matrix)
            ch.data[key] = elem[1:ch.iters, :]
        end
    end
end

"""
    nprods(ch)

Number of matrix-vector products plus transposed matrix-vector products
logged in `ConvergenceHistory` `ch`.
"""
nprods(ch::ConvergenceHistory) = ch.mvps+ch.mtvps

"""
    niters(ch)

Number of iterations logged in `ConvergenceHistory` `ch`.
"""
niters(ch::ConvergenceHistory) = ch.iters

"""
    nrests(ch)

Number of restarts logged in `ConvergenceHistory` `ch`.
"""
nrests(ch::RestartedHistory) = Int(ceil(ch.iters/ch.restart))

"""
    tolkeys(ch)

Key iterator of the tolerances logged in `ConvergenceHistory` `ch`.
"""
tolkeys(ch::ConvergenceHistory) = keys(ch.tol)

"""
    datakeys(ch)

Key iterator of the per iteration data logged in `ConvergenceHistory` `ch`.
"""
datakeys(ch::ConvergenceHistory) = keys(ch.data)

"""
    setmvps(ml, val)

Set `val` as number of matrix-vector products in [`MethodLog`](@ref) `ch`.
"""
setmvps(ch::ConvergenceHistory, val::Int) = ch.mvps=val

"""
    setmtvps(ml, val)

Set `val` as number of transposed matrix-vector products in [`MethodLog`](@ref) `ch`.
"""
setmtvps(ch::ConvergenceHistory, val::Int) = ch.mtvps=val

"""
    setconv(ml, val)

Set `val` as convergence status of the method in [`MethodLog`](@ref) ch.
"""
setconv(ch::ConvergenceHistory, val::Bool) = ch.isconverged=val

"""
    nextiter!(ml)

Adds one the the number of iterations in [`MethodLog`](@ref) `ch`. This is
necessary to avoid overwriting information with `push!(ml)`.
"""
nextiter!(ch::ConvergenceHistory) = ch.iters+=1

#### Errors
export PosSemidefException

type PosSemidefException <: Exception
    msg :: AbstractString
    PosSemidefException(msg::AbstractString="Matrix was not positive semidefinite") = new(msg)
end

#### Support for linear-operators-defined-as-functions
#
# Suppose you have a function implementing multiplication of b by A.
# This function can have any syntax, but for the purposes of
# illustration let's suppose it's defined as
#   mulbyA!(output, b, Adata)
# Where Adata might be some parameters that your function needs.
# You can represent it as a linear operator using
#   A = MatrixFcn{T}(m, n, (output, b) -> mulbyA!(output, b, Adata))
# where T is the "element type" of Adata.
# Note that there are a couple of requirements:
#   - mulbyA! stores the result in the pre-allocated output
#   - mulbyA! should also return output as its sole return value

# If the algorithm also needs multiplication by A', use MatrixCFcn
# instead and initialize it with both functions.

export MatrixFcn, MatrixCFcn

abstract AbstractMatrixFcn{T}

type MatrixFcn{T} <: AbstractMatrixFcn{T}
    m::Int
    n::Int
    mul::Function
end

type MatrixCFcn{T} <: AbstractMatrixFcn{T}
    m::Int
    n::Int
    mul::Function
    mulc::Function
end

MatrixFcn(A::AbstractMatrix) = MatrixFcn{eltype(A)}(size(A, 1), size(A, 2), (output, x)-> A_mul_B!(output, A, x))
MatrixCFcn(A::AbstractMatrix) = MatrixFcn{eltype(A)}(size(A, 1), size(A, 2), (output, x)->A_mul_B!(output, A, x), (output, b) -> Ac_mul_B(output, A, b))

eltype{T}(op::AbstractMatrixFcn{T}) = T

ndims(op::AbstractMatrixFcn) = 2

size(op::AbstractMatrixFcn, dim::Integer) = (dim == 1) ? op.m :
                                            (dim == 2) ? op.n : 1
size(op::AbstractMatrixFcn) = (op.m, op.n)

length(op::AbstractMatrixFcn) = op.m*op.n

*(op::AbstractMatrixFcn, b) = A_mul_B(op, b)

A_mul_B{R,S}(op::AbstractMatrixFcn{R}, b::AbstractVector{S}) = A_mul_B!(Array(promote_type(R,S), op.m), op, b)

A_mul_B!(output, op::AbstractMatrixFcn, b) = op.mul(output, b)

Ac_mul_B{R,S}(op::MatrixCFcn{R}, b::AbstractVector{S}) = Ac_mul_B!(Array(promote_type(R,S), op.n), op, b)

Ac_mul_B!(output, op::AbstractMatrixFcn, b) = op.mulc(output, b)
