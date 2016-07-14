using UnicodePlots
import  Base: last, keys, setindex!, getindex, \, eltype, empty!, eps, length,
        ndims, push!, real, size, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export A_mul_B, iters, last, Master

\(f::Function, b::Vector) = f(b)

#### Reporting
abstract MethodLog
type DummyHistory <: MethodLog end
type ConvergenceHistory{T} <: MethodLog
    mvps::Int
    mtvps::Int
    iters::Int
    restart::T
    isconverged::Bool
    tol::Dict{Symbol, Real}
    data::Dict{Symbol, Vector}
end
typealias PlainHistory ConvergenceHistory{Void}
typealias RestartedHistory ConvergenceHistory{Int}

DummyHistory() = DummyHistory()
ConvergenceHistory() = ConvergenceHistory(0,0,0,nothing,false,
                        Dict{Symbol, Real}(), Dict{Symbol, Vector}()
                        )
ConvergenceHistory(restart::Int) = ConvergenceHistory(0,0,0,restart,false,
                                    Dict{Symbol, Real}(), Dict{Symbol, Vector}()
                                    )

# Consulting

function getindex(ch::ConvergenceHistory, s::Symbol)
    haskey(ch.tol, s) && return ch.tol[s]
    ch.data[s]
end
getindex(ch::ConvergenceHistory, s::Symbol, i::Int) = ch.data[s][i]

products(ch::ConvergenceHistory) = ch.mvps+ch.mtvps

iters(ch::PlainHistory) = ch.iters
iters(ch::RestartedHistory) = ceil(ch.iters/ch.restart)

last(ch::ConvergenceHistory, key::Symbol) = ch.data[key][ch.iters]

tolkeys(ch::ConvergenceHistory) = keys(ch.tol)
datakeys(ch::ConvergenceHistory) = keys(ch.data)

# State change

setmvps(ch::DummyHistory, val::Int) = nothing
setmvps(ch::ConvergenceHistory, val::Int) = ch.mvps=val

setmtvps(ch::DummyHistory, val::Int) = nothing
setmtvps(ch::ConvergenceHistory, val::Int) = ch.mtvps=val

setconv(ch::DummyHistory, val::Bool) = nothing
setconv(ch::ConvergenceHistory, val::Bool) = ch.isconverged=val

function reserve!(ch::PlainHistory, s::Symbol, maxiter::Int; T::Type=Float64)
    ch.data[s] = Vector{T}(maxiter)
end
function reserve!(ch::RestartedHistory, s::Symbol, maxiter::Int; T::Type=Float64)
    ch.data[s] = Vector{T}(maxiter*ch.restart)
end
function reserve!(ch::PlainHistory, s::Symbol, maxiter::Int, size::Int; T::Type=Float64)
    aux = ch.data[s] = Vector{Vector{T}}(maxiter)
    for i in 1:length(aux)
        aux[i] = Vector{T}(size)
    end
end
function reserve!(ch::RestartedHistory, s::Symbol, maxiter::Int, size::Int; T::Type=Float64)
    aux = ch.data[s] = Vector{Vector{T}}(maxiter*ch.restart)
    for i in 1:length(aux)
        aux[i] = Vector{T}(size)
    end
end

setindex!(ch::DummyHistory, tol::Real, s::Symbol) = nothing
setindex!(ch::ConvergenceHistory, tol::Real, s::Symbol) = ch.tol[s] = tol

setindex!(ch::DummyHistory, val, s::Symbol, i::Int) = nothing
setindex!(ch::ConvergenceHistory, val, s::Symbol, i::Int) = ch.data[s][i] = val

setindex!(ch::DummyHistory, vec::Vector, s::Symbol) = nothing
setindex!(ch::ConvergenceHistory, vec::Vector, s::Symbol) = ch.data[s] = vec

nextiter!(::DummyHistory) = nothing
nextiter!(ch::ConvergenceHistory) = ch.iters+=1

push!(::DummyHistory, ::Symbol, ::Any) = nothing
push!(ch::ConvergenceHistory, key::Symbol, val) = ch.data[key][ch.iters] = val

shrink!(::DummyHistory) = nothing
shrink!(ch::ConvergenceHistory, key::Symbol) = resize!(ch.data[key], ch.iters)
function shrink!(ch::ConvergenceHistory)
    for key in datakeys(ch)
        shrink!(ch,key)
    end
end

#Plotting
function _plot{T<:Real}(vals::Vector{T}, iters::Int, gap::Int;
    restarts=ceil(iters/gap), color::Symbol=:blue, name::AbstractString="",
    title::AbstractString="", left::Int=1
    )
    maxy = maximum(vals)
    miny = minimum(vals)
    plot = lineplot([left],[miny],xlim=[left,iters],ylim=[miny,maxy],title=title,name=name)
    right = min(left+gap-1,iters)
    lineplot!(plot,collect(left:right),vals[left:right],color=color)
    for restart in 2:restarts
        left+=gap
        right = min(left+gap-1,iters)
        lineplot!(plot,collect(left:right),vals[left:right],color=color)
        lineplot!(plot,[left,left],[miny,maxy], color=:white)
    end
    plot
end

function _plot{T<:Real}(vals::Vector{Vector{T}}, iters::Int, gap::Int;
    restarts=ceil(iters/gap), color::Symbol=:blue, name::AbstractString="",
    title::AbstractString="", left::Int=1
    )
    n = size(vals[1],1)
    maxy = maximum(map(x -> maximum(x), vals))
    miny = minimum(map(x -> minimum(x), vals))
    plot = scatterplot([left],[miny],xlim=[left,iters],ylim=[miny,maxy],title=title,name=name)
    for i in left:iters
        scatterplot!(plot,[i for j in 1:n],vals[i][1:n],color=color)
    end
    for restart in 2:restarts
        left+=gap
        scatterplot!(plot,[left,left],[miny,maxy], color=:white)
    end
    plot
end

function showplot(ch::ConvergenceHistory; kwargs...)
    println("\n")
    for key in datakeys(ch)
        try
            restart = isa(ch, PlainHistory) ? ch.iters : ch.restart
            draw = _plot(ch.data[key], ch.iters, restart; name=string(key), kwargs...)
            println("$draw\n\n")
        catch e
            error("$e")
        end
    end
end

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

#Master type
typealias Master Tuple{Any,ConvergenceHistory}

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
