using UnicodePlots
import  Base: last, keys, setindex!, getindex, \, eltype, empty!, eps, length,
        ndims, push!, real, size, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export A_mul_B, iters, last, Master

\(f::Function, b::Vector) = f(b)

#Logging
abstract IterLog{T}
type SingleValue{T} <: IterLog{T}
    data::T
end
type MultiValue{T} <: IterLog{T}
    data::Vector{T}
end
SingleValue(T::Type) = SingleValue(zero(T))
MultiValue(T::Type, n::Int) = MultiValue(zeros(T,n))
getindex(il::IterLog, i::Int) = il.data[i]
setindex!(il::MultiValue, i::Int, val) = (il.data[i] = val)

#typealias ElemLog Vector{IterLog}
abstract MethodLog
type DummyLog <: MethodLog
    data::Dict{Symbol,IterLog}
end
MethodLog() = DummyLog(Dict{Symbol,IterLog}())

getindex(ml::DummyLog, key::Symbol) = ml.data[key]
setindex!(ml::DummyLog, val::IterLog, key::Symbol) = (ml.data[key]=val)

type MasterLog <: MethodLog
    iter::Int
    maxiter::Int
    restart::Int
    data::Dict{Symbol,Vector{IterLog}}
end
MethodLog(maxiter::Int) =
    MasterLog(0,maxiter,maxiter,Dict{Symbol, Vector{IterLog}}())
MethodLog(maxiter::Int, restart::Int) =
    MasterLog(0,maxiter*restart,restart,Dict{Symbol, Vector{IterLog}}())

getindex(ml::MasterLog, key::Symbol) = ml.data[key]
setindex!(ml::MasterLog, val::Vector{IterLog}, key::Symbol) = (ml.data[key]=val)
getindex(ml::MasterLog, key::Symbol, i::Int) = ml.data[key][i]
setindex!(ml::MasterLog, val::IterLog, key::Symbol, i::Int) = (ml.data[key][i]=val)

function add!(ml::MasterLog, key::Symbol; n::Int=1, T::Type=Float64)
    if n == 1
        ml.data[key] = Vector{SingleValue{T}}([SingleValue(T) for i in 1:ml.maxiter])
    else
        ml.data[key] = Vector{MultiValue{T}}([MultiValue(T,n) for i in 1:ml.maxiter])
    end
    nothing
end

iters(ml::MasterLog) = ml.iter

last(ml::DummyLog, key::Symbol) = ml[key].data
last(ml::MasterLog, key::Symbol) = ml[key][ml.iter].data

next!(ml::DummyLog) = nothing
next!(ml::MasterLog) = (ml.iter += 1; nothing)

push!(ml::DummyLog, key::Symbol, val) = (ml[key] = SingleValue(val))
function push!(ml::MasterLog, key::Symbol, val)
    ml[key][ml.iter].data = val
end

keys(ml::MethodLog) = keys(ml.data)

shrink!(ml::MasterLog, key::Symbol) = (resize!(ml[key], ml.iter); nothing)
function shrink!(ml::MasterLog)
    for key in keys(ml.data)
        shrink!(ml,key)
    end
end

isconverged(ml::MethodLog, key::Symbol, tol) = 0 <= last(ml, key) < tol

#Plotting
function plot{T<:Real}(vals::Vector{T}, iters::Int, gap::Int;
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

function plot{T<:Real}(vals::Vector{Vector{T}}, iters::Int, gap::Int;
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

showplot{T<:Real}(::Type{SingleValue{T}}, els::Vector{IterLog}, iters::Int, gap::Int; kwargs...) =
    println(plot(Real[els[i].data for i in 1:iters], iters, gap; kwargs...))

showplot{T<:Real}(::Type{MultiValue{T}}, els::Vector{IterLog}, iters::Int, gap::Int; kwargs...) =
    println(plot(Vector{Real}[els[i].data for i in 1:iters], iters, gap; kwargs...))

showplot(els::Vector{IterLog}, iters::Int, gap::Int; kwargs...) =
    showplot(typeof(els[1]), els::Vector{IterLog}, iters::Int, gap::Int; kwargs...)

function showplot(ml::MasterLog; kwargs...)
    println("\n")
    for key in keys(ml)
        try
            showplot(ml[key], ml.iter, ml.restart; name=string(key), kwargs...)
            println("\n\n")
        catch e
            warn("When trying to plot got the following error: $e")
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

#### Reporting
type ConvergenceHistory{T}
    isconverged::Bool
    threshold::T
    mvps::Int
    data::MethodLog
end

iters(ch::ConvergenceHistory) = iters(ch.data)

last(ch::ConvergenceHistory, key::Symbol) = last(ch.data,key)

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
