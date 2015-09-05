import Distributions: Categorical
export rkaczmarz, rkaczmarz!

##############################################################################
##
## "A randomized Kaczmarz algorithm with exponential convergence"
## Thomas Strohmer and Roman Vershynin
## At must be the transpose of A in Ax = b
## 
##############################################################################

function rkaczmarz(At, b; tol::Real=1e-10, maxiter::Integer = 100_000)
    T = Adivtype(At, b)
    x = zeros(T, size(At, 1))
    rkaczmarz!(x, At, b; tol = tol, maxiter = maxiter)
end

function rkaczmarz!(x, At, b; tol::Real=1e-10, maxiter::Integer = 100_000)

    # checks
    length(x) == size(At, 1) || error("The number of rows in At is $(size(At, 1)) but the length of x is $(length(x))")
    length(b) == size(At, 2) || error("The number of columns in At is $(size(At, 2)) but the length of b is $(length(b))")

    # initialize
    Tx = eltype(x)
    ch = ConvergenceHistory(false, tol, 0, Tx[])

    # precompute norm (distribution) and invnorm = 1/norm[i] 
    len_b = length(b)
    norm2 = _norm2cols(At)
    invnorm2 = 1 ./ norm2
    p = scale!(norm2, 1/sum(norm2))
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

        # check convergence every length(b)
        ssr = norm(res)
        push!(ch, ssr)
        if ssr < tol
            ch.isconverged = true
        end
    end
    return x, ch
end


##############################################################################
##
## Sparse Matrix
##
##############################################################################

if VERSION < v"0.4.0-dev+1307" 
    rowvals(S::Base.SparseMatrixCSC) = S.rowval
    nzrange(S::Base.SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)
end

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


