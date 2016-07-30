import  Base: last, keys, setindex!, getindex, \, eltype, empty!, eps, length,
        ndims, push!, real, size, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export  A_mul_B, niters, products, tolkeys, datakeys, nrestarts, setindex!,
        getindex, last, Master

\(f::Function, b::VecOrMat) = f(b)
*(f::Function, b::VecOrMat) = f(b)

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
    data::Dict{Symbol, VecOrMat}
end
typealias PlainHistory ConvergenceHistory{Void}
typealias RestartedHistory ConvergenceHistory{Int}

ConvergenceHistory() = ConvergenceHistory(0,0,0,nothing,false,
                        Dict{Symbol, Real}(), Dict{Symbol, VecOrMat}()
                        )
ConvergenceHistory(restart::Int) = ConvergenceHistory(0,0,0,restart,false,
                                    Dict{Symbol, Real}(), Dict{Symbol, VecOrMat}()
                                    )

# Consulting

function getindex(ch::ConvergenceHistory, s::Symbol)
    haskey(ch.tol, s) && return ch.tol[s]
    ch.data[s]
end
getindex(ch::ConvergenceHistory, s::Symbol, kwargs...) = ch.data[s][kwargs...]

products(ch::ConvergenceHistory) = ch.mvps+ch.mtvps

niters(ch::ConvergenceHistory) = ch.iters
nrestarts(ch::RestartedHistory) = Int(ceil(ch.iters/ch.restart))

tolkeys(ch::ConvergenceHistory) = keys(ch.tol)
datakeys(ch::ConvergenceHistory) = keys(ch.data)

# State change

setmvps(ch::DummyHistory, val::Int) = nothing
setmvps(ch::ConvergenceHistory, val::Int) = ch.mvps=val

setmtvps(ch::DummyHistory, val::Int) = nothing
setmtvps(ch::ConvergenceHistory, val::Int) = ch.mtvps=val

setconv(ch::DummyHistory, val::Bool) = nothing
setconv(ch::ConvergenceHistory, val::Bool) = ch.isconverged=val

function reserve!(ch::ConvergenceHistory, s::Symbol, maxiter::Int; T::Type=Float64)
    ch.data[s] = Vector{T}(maxiter)
end
function reserve!(ch::ConvergenceHistory, s::Symbol, maxiter::Int, size::Int; T::Type=Float64)
    ch.data[s] = Matrix{T}(maxiter,size)
end

setindex!(ch::DummyHistory, tol::Real, s::Symbol) = nothing
setindex!(ch::ConvergenceHistory, tol::Real, s::Symbol) = ch.tol[s] = tol

setindex!(ch::DummyHistory, val, s::Symbol, kwargs...) = nothing
setindex!(ch::ConvergenceHistory, val, s::Symbol, kwargs...) = ch.data[s][kwargs...] = val

nextiter!(::DummyHistory) = nothing
nextiter!(ch::ConvergenceHistory) = ch.iters+=1

push!(::DummyHistory, ::Symbol, ::Any) = nothing
function push!(ch::ConvergenceHistory, key::Symbol, vec::Union{Vector,Tuple})
    matrix = ch.data[key]
    width = size(matrix,2)
    base = (ch.iters-1)*width
    for i in 1:min(width,length(vec))
        matrix[base+i] = vec[i]
    end
end
push!(ch::ConvergenceHistory, key::Symbol, val) = ch.data[key][ch.iters] = val

shrink!(::DummyHistory) = nothing
function shrink!(ch::ConvergenceHistory)
    for key in datakeys(ch)
        elem = ch.data[key]
        if isa(elem, Vector)
            resize!(elem, ch.iters)
        elseif isa(elem, Matrix)
            ch.data[key] = elem[1:ch.iters, :]
        end
    end
end

#Recipes
@recipe function chef(ch::ConvergenceHistory)
    frame = 1
    n = count_plotables(ch)
    n > 0 || error("No plotables")
    layout := (count_plotables(ch), 1)
    for (name, draw) in ch.data
        plotable(draw) || continue
        @series begin
            isa(draw, Vector) && (seriestype := :line; label:="$name")
            isa(draw, Matrix) && (seriestype := :scatter; title:="$name"; label:="")
            subplot := frame
            draw
        end
        frame+=1
    end
end

@recipe function chef(ch::ConvergenceHistory, name::Symbol)
    draw = ch[name]
    plotable(draw) || error("Not plotable")
    isa(draw, Vector) && (seriestype := :line; label:="$name")
    isa(draw, Matrix) && (seriestype := :scatter; title:="$name"; label:="")
    @series begin
        draw
    end
end

function count_plotables(ch::ConvergenceHistory)
    candidates = collect(values(ch.data))
    length(filter(identity, map(plotable, candidates)))
end

plotable{T<:Real}(::VecOrMat{T}) = true
plotable(::Any) = false

#Plotting
function _plot{T<:Real}(vals::Vector{T}, iters::Int, gap::Int;
    restarts=Int(ceil(iters/gap)), color::Symbol=:blue, name::AbstractString="",
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

function _plot{T<:Real}(vals::Matrix{T}, iters::Int, gap::Int;
    restarts=Int(ceil(iters/gap)), color::Symbol=:blue, name::AbstractString="",
    title::AbstractString="", left::Int=1
    )
    n = size(vals,2)
    maxy = maximum(vals)
    miny = minimum(vals)
    plot = scatterplot([left],[miny],xlim=[left,iters],ylim=[miny,maxy],title=title,name=name)
    for i in left:iters
        scatterplot!(plot,[i for j in 1:n],vals[i,1:n],color=color)
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
        catch
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
