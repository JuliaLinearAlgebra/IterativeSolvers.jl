import Distributions: Categorical
export rkaczmarz, rkaczmarz!

##############################################################################
##
## "A randomized Kaczmarz algorithm with exponential convergence"
## Thomas Strohmer and Roman Vershynin
## rkaczmarzT takes At as first argument to solve Ax = b
##############################################################################

function rkaczmarz(A, b ; kwargs...)
    rkaczmarzT(A', b ; kwargs...)
end

function rkaczmarzT(At, b; kwargs...)
    T = Adivtype(At, b)
    x = zeros(T, size(At, 1))
    rkaczmarzT!(x, At, b; kwargs...)
end

function rkaczmarzT!(x, At, b; tol::Real=1e-10, maxiter::Integer = 100_000)
    # checks
    length(x) == size(At, 1) || error("The number of rows in At is $(size(At, 1)) but the length of x is $(length(x))")
    length(b) == size(At, 2) || error("The number of columns in At is $(size(At, 2)) but the length of b is $(length(b))")

    # initialize
    Tx = eltype(x)
    ssrs = Tx[]

    # precompute norm (distribution) and invnorm = 1/norm[i] 
    len_b = length(b)
    norm2 = _norm2cols(At)
    invnorm2 = 1 ./ norm2
    @show any(x -> x == 0, norm2)
    any(x -> x == 0, norm2) && error("Matrix is singular")
    p = scale!(norm2, inv(sum(norm2)))
   	distribution = Categorical(p)
    res = b - At_mul_B(At, x)
    iter = 0
    while iter < maxiter 
        iter += 1
        for _ in 1:len_b                
            # draw a row
            i = rand(distribution)
            resi = b[i] - _dot(x, At, i)
            res[i] = resi
            α = resi * invnorm2[i]
            _update!(x, At, i, α)
        end

        # check convergence every length(b) iteration
        ssr = norm(res)
        push!(ssrs, ssr)
        if ssr < tol^2
            break
        end
    end
    return x, ConvergenceHistory(ssrs[end] < tol, tol, length(ssrs), ssrs)
end


##############################################################################
##
## Sparse Matrix
##
##############################################################################
function _norm2cols{Tv, Ti}(At::SparseMatrixCSC{Tv, Ti})
    Tv[sumabs2(sub(nonzeros(At),  nzrange(At, i))) for i in 1:size(At, 2)]
end

function _dot(x, At::SparseMatrixCSC, i::Integer)
    Atrows = rowvals(At)
    Atvals = nonzeros(At)
    T = Amultype(At, x)
    out = zero(T)
    @inbounds for k in nzrange(At, i)
        out += x[Atrows[k]] * Atvals[k]
    end
    return out
end

function _update!(x, At::SparseMatrixCSC, i::Integer, α)
    Atrows = rowvals(At)
    Atvals = nonzeros(At)
    @inbounds for k in nzrange(At, i)
        x[Atrows[k]] += α * Atvals[k]
    end
end

##############################################################################
##
## Full Matrix
##
##############################################################################
_norm2cols(At::Matrix) = vec(sumabs2(At, 1))

_dot(x, At::Matrix, i::Integer) = dot(x, slice(At, :, i))

function _update!(x, At::Matrix, i::Integer, α)
    @inbounds for k in 1:size(At, 1)
        x[k] += α * At[k, i]
    end
end


