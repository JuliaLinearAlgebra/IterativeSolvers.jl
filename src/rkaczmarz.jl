import Distributions: Categorical
export rkaczmarz, rkaczmarz!

##############################################################################
##
## "A randomized Kaczmarz algorithm with exponential convergence"
## Thomas Strohmer and Roman Vershynin
##
##############################################################################

function rkaczmarz(At, b; tol::Real=1e-10, maxiter::Integer = 100_000)
    T = Adivtype(At, b)
    x = zeros(T, size(At, 1))
    converged = rkaczmarz!(x, At, b; tol = tol, maxiter = maxiter)
    return x, converged
end

function rkaczmarz!(x, At, b; tol::Real=1e-10, maxiter::Integer = 100_000)

    length(x) == size(At, 1) || error("The number of rows in At is $(size(At, 1)) but the length of x is $(length(x))")
    length(b) == size(At, 2) || error("The number of columns in At is $(size(At, 2)) but the length of b is $(length(b))")

    # precompute norm (distribution) and invnorm = 1/norm[i] 
    len_b = length(b)
    norm2 = _norm2cols(At)
    invnorm2 = 1 ./ norm2
    p = scale!(norm2, 1/sum(norm2))
   	distribution = Categorical(p)

    iter = 0
    while iter < maxiter 
        iter += 1
        for _ in 1:len_b                
           
            # draw a row
            i = rand(distribution)
       
            α = (b[i] - _dot(x, At, i)) * invnorm2[i]

            # update x_k
            _update!(x, At, i, α)
        end

        # check maxabs(Ax-b) 
        if maxabs(x, At, b, tol)
            return true
        end
    end
    return false
end

function maxabs(x, At, b, tol::Real)
    for i in 1:length(b)
        current = b[i] - _dot(x, At, i)
        if abs(current) > tol
            return false
        end
    end
    return true
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

function _norm2cols(At::SparseMatrixCSC)
    Tv[sumabs2(sub(nonzeros(At),  nzrange(At, i))) for i in 1:size(At, 2)]
end

function _dot(x, At::SparseMatrixCSC, i::Integer)
    Atrows = rowvals(At)
    Atvals = nonzeros(At)
    T = Amultype(At, x)
    out = zero(T)
    for k in nzrange(At, i)
        out += x[Atrows[k]] * Atvals[k]
    end
    return out
end

function _update!(x, At::SparseMatrixCSC, i::Integer, α)
    Atrows = rowvals(At)
    Atvals = nonzeros(At)
    for k in nzrange(At, i)
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
    for k in 1:size(At, 1)
        x[k] += α * At[k, i]
    end
end


