import  Base: last, keys, setindex!, getindex, \, eltype, empty!, eps, length,
        ndims, push!, real, size, ctranspose, *, A_mul_B!, Ac_mul_B, Ac_mul_B!
export  A_mul_B, niters, nprods, tolkeys, datakeys, nrests, setindex!,
        getindex, last

# Improve readability of iterative methods
\(f::Function, b::VecOrMat) = f(b)
*(f::Function, b::VecOrMat) = f(b)

#### Reporting
"""
Log an iterative method's general and per-iteration information.
"""
abstract MethodLog

"""
Placeholder to put inside an iterative method instead of `ConvergenceHistory`
for not wasting memory. Doesn't actually store any kind of information.

**Implements**

* `Base`: `setindex!`, `push!`

"""
type DummyHistory <: MethodLog end

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
type ConvergenceHistory{T} <: MethodLog
    mvps::Int
    mtvps::Int
    iters::Int
    restart::T
    isconverged::Bool
    tol::Dict{Symbol, Real}
    data::Dict{Symbol, VecOrMat}
end
ConvergenceHistory() = ConvergenceHistory(0,0,0,nothing,false,
                        Dict{Symbol, Real}(), Dict{Symbol, VecOrMat}()
                        )
ConvergenceHistory(restart::Int) = ConvergenceHistory(0,0,0,restart,false,
                                    Dict{Symbol, Real}(), Dict{Symbol, VecOrMat}()
                                    )

"""
`ConvergeHistory` without resets.
"""
typealias PlainHistory ConvergenceHistory{Void}

"""
`ConvergeHistory` with resets.
"""
typealias RestartedHistory ConvergenceHistory{Int}

function getindex(ch::ConvergenceHistory, s::Symbol)
    haskey(ch.tol, s) && return ch.tol[s]
    ch.data[s]
end
getindex(ch::ConvergenceHistory, s::Symbol, kwargs...) = ch.data[s][kwargs...]

setindex!(ch::DummyHistory, tol::Real, s::Symbol) = nothing
setindex!(ch::ConvergenceHistory, tol::Real, s::Symbol) = ch.tol[s] = tol

setindex!(ch::DummyHistory, val, s::Symbol, kwargs...) = nothing
setindex!(ch::ConvergenceHistory, val, s::Symbol, kwargs...) = ch.data[s][kwargs...] = val

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

"""
    reserve!(ch, key, maxiter)
    reserve!(ch, key, maxiter, size)

Reserve space for per iteration data in `ch`. If size is provided, intead of a
vector it will reserve matrix of dimensions `(maxiter, size)`.

**Arguments**

* `ch::ConvergenceHistory`: convergence history.

* `key::Symbol`: key used to identify the data.

* `maxiter::Int`: number of iterations to save space for.

* `size::Int`: number of elements to store with the `key` identifier.

*Keywords*

* `T::Type=Float64`: type of the elements to store.

"""
function reserve!(ch::ConvergenceHistory, key::Symbol, maxiter::Int; T::Type=Float64)
    ch.data[key] = Vector{T}(maxiter)
end
function reserve!(ch::ConvergenceHistory, key::Symbol, maxiter::Int, size::Int; T::Type=Float64)
    ch.data[key] = Matrix{T}(maxiter,size)
end

"""
    shrink!(ml)

shrinks the reserved space for `MethodLog` `ch` to the space actually used to log.
"""
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
setmvps(ch::DummyHistory, val::Int) = nothing
setmvps(ch::ConvergenceHistory, val::Int) = ch.mvps=val

"""
    setmtvps(ml, val)

Set `val` as number of transposed matrix-vector products in [`MethodLog`](@ref) `ch`.
"""
setmtvps(ch::DummyHistory, val::Int) = nothing
setmtvps(ch::ConvergenceHistory, val::Int) = ch.mtvps=val

"""
    setconv(ml, val)

Set `val` as convergence status of the method in [`MethodLog`](@ref) ch.
"""
setconv(ch::DummyHistory, val::Bool) = nothing
setconv(ch::ConvergenceHistory, val::Bool) = ch.isconverged=val

"""
    nextiter!(ml)

Adds one the the number of iterations in [`MethodLog`](@ref) `ch`. This is
necessary to avoid overwriting information with `push!(ml)`.
"""
nextiter!(::DummyHistory) = nothing
nextiter!(ch::ConvergenceHistory) = ch.iters+=1


#Recipes
#See Plots.jl tutorial on recipes
#Plot entire ConvergenceHistory
@recipe function chef(ch::ConvergenceHistory; sep = :white)
    candidates = collect(values(ch.data))
    plotables = map(plotable, candidates)
    n = length(filter(identity, plotables))
    n > 0 || error("No plotable")
    frame = 1
    layout := (n, 1)
    for (name, draw) in collect(ch.data)[plotables]
        @series begin
            isa(draw, Vector) && (seriestype:= :line; label:="$name")
            isa(draw, Matrix) && (seriestype:= :scatter; title:="$name"; label:="")
            subplot := frame
            map(x->isnan(x) ? typeof(x)(0) : x,draw)
        end
        if isa(ch, RestartedHistory)
            label := ""
            linecolor := sep

            left=1
            maxy = round(maximum(draw),2)
            miny = round(minimum(draw),2)
            for restart in 2:nrests(ch)
                @series begin
                    left+=ch.restart
                    subplot := frame
                    [left,left],[miny,maxy]
                end
            end
        end
        frame+=1
    end
end

#Plot collection `ch[name]`.
@recipe function chef(ch::ConvergenceHistory, name::Symbol; sep = :white)
    draw = ch[name]
    plotable(draw) || error("Not plotable")
    isa(draw, Vector) && (seriestype-->:line; label-->"$name")
    isa(draw, Matrix) && (seriestype-->:scatter; title-->"$name"; label-->"")
    @series begin
        draw
    end
    if isa(ch, RestartedHistory)
        label := ""
        linecolor := sep

        left=1
        maxy = round(maximum(draw),2)
        miny = round(minimum(draw),2)
        for restart in 2:nrests(ch)
            @series begin
                left+=ch.restart
                [left,left],[miny,maxy]
            end
        end
    end
end

"""
    plotable(x)

Determine whether a collection `x` is plotable. Only vectors and matrices are
such objects.
"""
plotable{T<:Real}(::VecOrMat{T}) = true
plotable(::Any) = false

"""
    showplot(ch)

Print all plotable information inside `ConvergenceHistory` `ch`.

*Note:* This is what is called when the `plot` keyword is set.
"""
function showplot(ch::ConvergenceHistory)
    candidates = collect(values(ch.data))
    plotables = map(plotable, candidates)
    n = length(filter(identity, plotables))
    n > 0 || return
    println("\n")
    for (name, draw) in collect(ch.data)[plotables]
        restart = isa(ch, PlainHistory) ? ch.iters : ch.restart
        drawing = plot_collection(draw, ch.iters, restart; name=string(name))
        println("$drawing\n\n")
    end
end

"""
    plot_collection(x)

Build a `UnicodePlot.Plot` object from the plotable collection `x`.

If `x` is a vector, a series will be made. In case of being a matrix an scatterplot
will be returned.
"""
function plot_collection{T<:Real}(vals::Vector{T}, iters::Int, gap::Int;
    restarts=Int(ceil(iters/gap)), color::Symbol=:blue, name::AbstractString="",
    title::AbstractString="", left::Int=1
    )
    maxy = round(maximum(vals),2)
    miny = round(minimum(vals),2)
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
function plot_collection{T<:Real}(vals::Matrix{T}, iters::Int, gap::Int;
    restarts=Int(ceil(iters/gap)), color::Symbol=:blue, name::AbstractString="",
    title::AbstractString="", left::Int=1
    )
    n = size(vals,2)
    maxy = round(maximum(vals),2)
    miny = round(minimum(vals),2)
    plot = scatterplot([left],[miny],xlim=[left,iters],ylim=[miny,maxy],title=title,name=name)
    for i in left:iters
        scatterplot!(plot,[i for j in 1:n],vec(vals[i,1:n]),color=color)
    end
    for restart in 2:restarts
        left+=gap
        lineplot!(plot,[left,left],[miny,maxy], color=:white)
    end
    plot
end

#### Type-handling
"""
    Adivtype(A, b)

Determine type of the division of an element of `b` against an element of `A`:

`typeof(one(eltype(b))/one(eltype(A)))`
"""
Adivtype(A, b) = typeof(one(eltype(b))/one(eltype(A)))

"""
    Amultype(A, x)

Determine type of the multiplication of an element of `b` with an element of `A`:

`typeof(one(eltype(A))*one(eltype(x)))`
"""
Amultype(A, x) = typeof(one(eltype(A))*one(eltype(x)))
if VERSION < v"0.4.0-dev+6068"
    real{T<:Real}(::Type{Complex{T}}) = T
    real{T<:Real}(::Type{T}) = T
end

eps{T<:Real}(::Type{Complex{T}}) = eps(T)

"""
    randx(A, b)

Build a random unitary vector `Vector{T}`, where `T` is `Adivtype(A,b)`.
"""
function randx(A, b)
    T = Adivtype(A, b)
    x = initrand!(Array(T, size(A, 2)))
end

"""
    zerox(A, b)

Build a zeros vector `Vector{T}`, where `T` is `Adivtype(A,b)`.
"""
function zerox(A, b)
    T = Adivtype(A, b)
    x = zeros(T, size(A, 2))
end

#### Numerics
"""
    update!(x, α::Number, p::AbstractVector)

Overwrite `x` with the result of `x += α*p`.
"""
function update!(x, α::Number, p::AbstractVector)
    for i = 1:length(x)
        x[i] += α*p[i]
    end
    x
end

"""
    initrand!(v)

Overwrite `v` with a random unitary vector of the same length.
"""
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

#### Errors
export PosSemidefException

type PosSemidefException <: Exception
    msg :: AbstractString
    PosSemidefException(msg::AbstractString="Matrix was not positive semidefinite") = new(msg)
end

export FuncMat

"""
    FuncMat{T}

**Fields**

* `m::Int` = number of columns.

* `n::Int` = number of rows.

* `mul::Function` = `A*b` implementation.

* `cmul::Function` = `A'*b` implementation.

**Constructors**

    FuncMat(A)
    FuncMat(
        m::Int, n::Int; typ::Type=Float64, ctrans::Bool=false,
        mul::Function=identity, cmul=identity
    )

**Arguments**

* `A::AbstractMatrix` = matrix.

* `m::Int` = number of columns.

* `n::Int` = number of rows.

* `typ::Type = Float64` = `eltype(::FuncMat)`.

* `mul::Function = identity` = `A*b` implementation.

* `cmul::Function = identity` = `A'*b` implementation.

"""
type FuncMat{T}
    m::Int
    n::Int
    mul::Function
    cmul::Function
end
function FuncMat(
        m::Int, n::Int; typ::Type=Float64,
        mul::Function=identity, cmul=mul
        )
    FuncMat{typ}(m, n, mul, cmul)
end
FuncMat(A::AbstractMatrix) = FuncMat(
    size(A, 1),
    size(A, 2),
    typ = eltype(A),
    mul = (output, x)->A_mul_B!(output, A, x),
    cmul = (output, b) -> Ac_mul_B(output, A, b)
    )

eltype{T}(::FuncMat{T}) = T

ndims(::FuncMat) = 2

size(op::FuncMat) = (op.m, op.n)
size(op::FuncMat, dim::Integer) = (dim == 1) ? op.m : (dim == 2) ? op.n : 1

length(op::FuncMat) = op.m*op.n

ctranspose{T}(op::FuncMat{T}) = FuncMat{T}(op.n, op.m, op.cmul, op.mul)

*(op::FuncMat, b) = A_mul_B(op, b)

function A_mul_B{R,S}(op::FuncMat{R}, b::AbstractVector{S})
    A_mul_B!(Array(promote_type(R,S), op.m), op, b)
end
function A_mul_B{R,S}(op::FuncMat{R}, b::AbstractMatrix{S})
    A_mul_B!(Array(promote_type(R,S), op.m, size(b,2)), op, b)
end

function A_mul_B!(output, op::FuncMat, b::AbstractVector)
    op.mul == identity && error("A*b not defined")
    op.mul(output, b)
end
function A_mul_B!(output, op::FuncMat, b::AbstractMatrix)
    op.mul == identity && error("A*b not defined")
    columns = [op.mul(output, b[:,i]) for i in 1:size(b,2)]
    hcat(columns...)
end

function Ac_mul_B{R,S}(op::FuncMat{R}, b::AbstractVector{S})
    Ac_mul_B!(Array(promote_type(R,S), op.n), op, b)
end
function Ac_mul_B{R,S}(op::FuncMat{R}, b::AbstractMatrix{S})
    Ac_mul_B!(Array(promote_type(R,S), op.n, size(b,2)), op, b)
end

function Ac_mul_B!(output, op::FuncMat, b::AbstractVector)
    op.cmul == identity && error("A'*b not defined")
    op.cmul(output, b)
end
function Ac_mul_B!(output, op::FuncMat, b::AbstractMatrix)
    op.mul == identity && error("A'*b not defined")
    columns = [op.cmul(output, b[:,i]) for i in 1:size(b,2)]
    hcat(columns...)
end
