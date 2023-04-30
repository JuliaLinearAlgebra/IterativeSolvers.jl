using RecipesBase

import Base: getindex, setindex!, push!, keys, show

export ConvergenceHistory
export nprods, niters, nrests

########
# Type #
########

"""
Store general and in-depth information about an iterative method.

# Fields

`mvps::Int`: number of matrix vector products.

`mtvps::Int`: number of transposed matrix-vector products

`iters::Int`: iterations taken by the method.

`restart::T`: restart relevant information.

* `T == Int`: iterations per restart.
* `T == Nothing`: methods without restarts.

`isconverged::Bool`: convergence of the method.

`data::Dict{Symbol,Any}`: Stores all the information stored during the method execution.
It stores tolerances, residuals and other information, e.g. Ritz values in [`svdl`](@ref).

# Constructors

    ConvergenceHistory()
    ConvergenceHistory(restart)

Create `ConvergenceHistory` with empty fields.

# Arguments

`restart`: number of iterations per restart.

# Plots

Supports plots using the `Plots.jl` package via a type recipe. Vectors are
ploted as series and matrices as scatterplots.

# Implements

`Base`: `getindex`, `setindex!`, `push!`

"""
mutable struct ConvergenceHistory{T,K}
    mvps::Int
    mtvps::Int
    iters::Int
    restart::K
    isconverged::Bool
    data::Dict{Symbol, Any}
end
function ConvergenceHistory(;restart=nothing, partial=false)
    ConvergenceHistory{!partial,typeof(restart)}(
        0,0,0,restart,false,Dict{Symbol, Any}()
        )
end

"""
Stores information of the current iteration.
"""
const PartialHistory = ConvergenceHistory{false}

"""
Stores the information of all the iterations.
"""
const CompleteHistory = ConvergenceHistory{true}

"""
History without resets.
"""
const UnrestartedHistory{T} = ConvergenceHistory{T, Nothing}

"""
History with resets.
"""
const RestartedHistory{T} = ConvergenceHistory{T, Int}


#############
# Functions #
#############

function show(io::IO, ch::ConvergenceHistory)
    print(io, ch.isconverged ? "Converged" : "Not converged",
        " after ", ch.iters, " iterations.")
end

"""
    getindex(ch, s)

Get collection or tolerance associated with key `s` in `ch::ConvergenceHistory`.

    getindex(ch, s, kwargs...)

Access elements of the collection associated with key `s` in `ch::ConvergenceHistory`.
"""
getindex(ch::ConvergenceHistory, s::Symbol) = ch.data[s]
getindex(ch::ConvergenceHistory, s::Symbol, kwargs...) = ch.data[s][kwargs...]

"""
    setindex!(ch, tol, s)

Set tolerance value associated with `s` in `ch::ConvergenceHistory` to `tol`.

    setindex!(ch, val, s, kwargs...)

Set collection element associated with key `s` in `ch::ConvergenceHistory` to val.
"""
setindex!(ch::ConvergenceHistory, val, s::Symbol) = ch.data[s] = val
setindex!(ch::ConvergenceHistory, val, s::Symbol, kwargs...) = ch.data[s][kwargs...] = val

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
push!(ch::ConvergenceHistory, key::Symbol, data) = push_custom_data!(ch, key, data)

push_custom_data!(ch::PartialHistory, key::Symbol, data) = ch.data[key] = data
push_custom_data!(ch::CompleteHistory, key::Symbol, data) = ch.data[key][ch.iters] = data

"""
    reserve!(ch, key, maxiter)
    reserve!(typ, key, maxiter)
    reserve!(ch, key, maxiter, size)
    reserve!(typ, ch, key, maxiter, size)

Reserve space for per iteration data in `ch`. If size is provided, instead of a
vector it will reserve matrix of dimensions `(maxiter, size)`.

# Arguments

`typ::Type`: Type of the elements to store. Defaults to `Float64` when not given.

`ch::ConvergenceHistory`: convergence history.

`key::Union{Symbol,Vector{Symbol}}`: key used to identify the data.

`maxiter::Int`: number of iterations to save space for.

`size::Int`: number of elements to store with the `key` identifier.

"""
function reserve!(ch::ConvergenceHistory, keys::Vector{Symbol}, kwargs...)
    for key in keys
        reserve!(ch, key, kwargs...)
    end
end
function reserve!(ch::ConvergenceHistory, key::Symbol, kwargs...)
    reserve!(Float64, ch, key, kwargs...)
end
function reserve!(typ::Type, ch::ConvergenceHistory, key::Symbol, kwargs...)
    _reserve!(typ, ch, key, kwargs...)
end

#If partialhistory, theres no need to store a vector or matrix, instead
#store nothing or store a vector respectively.
_reserve!(typ::Type, ch::PartialHistory, key::Symbol, ::Int) = nothing
function _reserve!(typ::Type, ch::PartialHistory, key::Symbol, ::Int, size::Int)
    ch.data[key] = Vector{typ}(undef, size)
end
function _reserve!(typ::Type, ch::CompleteHistory, key::Symbol, len::Int)
    ch.data[key] = Vector{typ}(undef, len)
end
function _reserve!(typ::Type, ch::CompleteHistory, key::Symbol, len::Int, size::Int)
    ch.data[key] = Matrix{typ}(undef, len, size)
end

"""
    shrink!(ml)

shrinks the reserved space for `ConvergenceHistory` `ch` to the space actually used to log.
"""
shrink!(::PartialHistory) = nothing
function shrink!(ch::CompleteHistory)
    for key in keys(ch)
        elem = ch.data[key]
        if isa(elem, Vector)
            resize!(elem, ch.iters)
        elseif isa(elem, Matrix)
            ch.data[key] = elem[1:ch.iters, :]
        end
    end
end

"""
    nextiter!(ml)

Adds one the the number of iterations in [`ConvergenceHistory`](@ref) `ch`. This is
necessary to avoid overwriting information with `push!(ml)`. It is also able
to update other information of the method.
"""
function nextiter!(ch::ConvergenceHistory; mvps=0,mtvps=0)
    ch.iters+=1
    ch.mvps+=mvps
    ch.mtvps+=mtvps
end

"""
    keys(ch)

Key iterator of the per iteration data logged in `ConvergenceHistory` `ch`.
"""
keys(ch::ConvergenceHistory) = keys(ch.data)

"""
    setconv(ml, val)

Set `val` as convergence status of the method in [`ConvergenceHistory`](@ref) ch.
"""
setconv(ch::ConvergenceHistory, val::Bool) = ch.isconverged=val

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

#########
# Plots #
#########

"""
    plotable(x)

Determine whether a collection `x` is plotable. Only vectors and matrices are
such objects.
"""
plotable(::VecOrMat{T}) where {T <: Real} = true
plotable(::Any) = false


# Plot entire ConvergenceHistory. `sep` is the color of the restart separator.
@recipe function chef(ch::CompleteHistory; sep = :white)
    candidates = collect(values(ch.data))
    plotables = convert(Vector{Bool}, map(plotable, candidates))
    n = length(filter(identity, plotables))
    n > 0 || error("No plotables")
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
            maxy = round(maximum(draw); digits=2)
            miny = round(minimum(draw); digits=2)
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

# Plot collection `ch[name]`. `sep` is the color of the restart separator.
@recipe function chef(ch::CompleteHistory, name::Symbol; sep = :white)
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
        maxy = round(maximum(draw); digits=2)
        miny = round(minimum(draw); digits=2)
        for restart in 2:nrests(ch)
            @series begin
                left+=ch.restart
                [left,left],[miny,maxy]
            end
        end
    end
end
