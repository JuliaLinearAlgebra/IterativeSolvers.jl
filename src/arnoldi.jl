export eigvals_arnoldi

## Arnoldi factorization starting from step-k
function arnoldifac!{T}(k::Int, m::Int, V::Matrix{T}, H::Matrix{T}, f::Vector{T}, A, prec::T)
    if m <= k
        error("m must be greater than k")
    end

    ## Keep the upperleft k x k submatrix of H and
    ## set other elements to 0
    H[:, (k + 1):end] = zero(T)
    H[(k + 1):end, 1:k] = zero(T)

    beta::T = norm(f)

    for i = k:(m - 1)
        ## V{i+1} <- f / ||f||
        V[:, i + 1] = f / beta
        H[i + 1, i] = beta

        ## w <- A * V{i+1}
        w::Vector{T} = A * V[:, i + 1]

        H[i, i + 1] = beta
        H[i + 1, i + 1] = dot(V[:, i + 1], w)

        ## f <- w - V * V' * w
        f[:] = w - beta * V[:, i] - H[i + 1, i + 1] * V[:, i + 1]
        beta = norm(f)

        ## f/||f|| is going to be the next column of V, so we need to test
        ## whether V' * (f/||f||) ~= 0
        Vf::Vector{T} = V[:, 1:(i + 1)]' * f
        count = 0
        while count < 5 && maximum(abs(Vf)) > prec * beta
            ## f <- f - V * Vf
            f[:] -= V[:, 1:(i + 1)] * Vf
            ## h <- h + Vf
            H[i, i + 1] += Vf[i]
            H[i + 1, i] = H[i, i + 1]
            H[i + 1, i + 1] += Vf[i + 1]
            ## beta <- ||f||
            beta = norm(f)

            Vf[:] = V[:, 1:(i + 1)]' * f
            count += 1
        end
    end
end

## Apply shifts on V, H and f
function applyshifts!{T}(k::Int, V::Matrix{T}, H::Matrix{T}, f::Vector{T}, shifts::Vector{T})
    n = size(V, 1)
    ncv = size(V, 2)
    Q::Matrix{T} = eye(T, ncv)

    for i = (k + 1):ncv
        ## QR decomposition of H-mu*I, mu is the shift
        for j = 1:ncv
            H[j, j] -= shifts[i]
        end
        qr = tridiagqr!(H)  ## H -> RQ
        applyright!(qr, Q)  ## Q -> Q * Qi
        for j = 1:ncv
            H[j, j] += shifts[i]
        end
    end

    ## V -> VQ, only need to update the first k+1 columns
    ## Q has some elements being zero
    ## The first (ncv - k + i) elements of the i-th column of Q are non-zero
    # Vs = zeros(T, n, k + 1)
    # for i = 1:k
    #     nnz = ncv - k + i
    #     Vs[:, i] = V[:, 1:nnz] * Q[1:nnz, i]
    # end
    # Vs[:, k + 1] = V * Q[:, k + 1]
    # V[:, 1:(k + 1)] = Vs
    ## However this seems to be slower than a direct multiplication

    ## V -> VQ
    V[:, :] = V * Q
    f[:] = f * Q[ncv, k] + V[:, k + 1] * H[k + 1, k]
end

## Retrieve Ritz values and Ritz vectors
function ritzpairs!{T}(which::Symbol, H::Matrix{T}, ritzval::Vector{T}, ritzvec::Matrix{T}, ritzest::Vector{T})
    ncv = size(ritzvec, 1)
    neigs = size(ritzvec, 2)
    ## Eigen decomposition on H, which is symmetric and tridiagonal
    decomp = eigfact(SymTridiagonal(diag(H), diag(H, -1)))

    ## Sort Ritz values according to "which"
    if which == :LM
        trans = abs
        rev = true
    elseif which == :SM
        trans = abs
        rev = false
    elseif which == :LA || which == :BE
        trans = (x -> x)
        rev = true
    else which == :SA
        trans = (x -> x)
        rev = false
    end

    ix = sortperm(decomp.values, by = trans, rev = rev)

    if which == :BE
        ixcp = copy(ix)
        for i = 1:ncv
            if i % 2 == 1
                ix[i] = ixcp[(i + 1) / 2]
            else
                ix[i] = ixcp[ncv - i / 2 + 1]
            end
        end
    end

    ritzval[:] = decomp.values[ix]
    ritzest[:] = decomp.vectors[ncv, ix]
    ritzvec[:, :] = decomp.vectors[:, ix[1:neigs]]
end

## Adjusted neigs
function neigsadjusted{T}(neigs::Int, ncv::Int, nconv::Int, ritzest::Vector{T}, prec::T)
    neigsnew = neigs + sum(abs(ritzest[(neigs + 1):ncv]) .< prec)
    neigsnew += min(nconv, div(ncv - neigsnew, 2))
    if neigsnew == 1 && ncv >= 6
        neigsnew = div(ncv, 2)
    elseif neigsnew == 1 && ncv > 2
        neigsnew = 2
    end

    return neigsnew
end

function eigvals_arnoldi(A, neigs::Int = 6;
                         ncv::Int = min(size(A, 1), max(2 * neigs + 1, 20)),
                         which::Symbol = :LM,
                         tol::Real = sqrt(eps(eltype(A))),
                         maxiter::Int = 300,
                         returnvec::Bool = true,
                         v0 = rand(eltype(A), size(A, 1)) - convert(eltype(A), 0.5))
    T = eltype(A)
    ## Size of matrix
    n = size(A, 1)
    ## Check arguments
    if neigs < 1 || neigs > n - 1
        error("neigs must satisfy 1 <= neigs <= n - 1, n is the size of matrix")
    end
    if ncv <= neigs || ncv > n
        error("ncv must satisfy neigs < ncv <= n, n is the size of matrix")
    end

    ## Number of matrix operations called
    nmatop::Int = 0
    ## Number of restarting iteration
    niter::Int = maxiter
    ## Number of converged eigenvalues
    nconv::Int = 0

    ## Matrices and vectors in the Arnoldi factorization
    V::Matrix{T} = zeros(T, n, ncv)
    H::Matrix{T} = zeros(T, ncv, ncv)
    f::Vector{T} = zeros(T, n)
    ritzval::Vector{T} = zeros(T, ncv)
    ritzvec::Matrix{T} = zeros(T, ncv, neigs)
    ritzest::Vector{T} = zeros(T, ncv)
    ritzconv::Vector{Bool} = zeros(Bool, neigs)

    ## Precision parameter used to test convergence
    prec::T = eps(T)^(convert(T, 2.0) / 3)

    ## Initialize vectors
    v0norm = norm(v0)
    if v0norm < prec
        error("initial residual vector cannot be zero")
    end
    v0[:] /= v0norm
    w::Vector{T} = A * v0
    nmatop += 1
    H[1, 1] = dot(v0, w)
    f[:] = w - v0 * H[1, 1]
    V[:, 1] = v0

    ## First Arnoldi factorization
    arnoldifac!(1, ncv, V, H, f, A, prec)
    nmatop += (ncv - 1)
    ritzpairs!(which, H, ritzval, ritzvec, ritzest)
    ## Restarting
    thresh::Vector{T} = zeros(T, neigs)
    resid::Vector{T} = zeros(T, neigs)
    for i = 1:maxiter
        ## Convergence test
        thresh[:] = tol * max(abs(ritzval[1:neigs]), prec)
        resid[:] = abs(ritzest[1:neigs]) * norm(f)
        ritzconv[:] = (resid .< thresh)
        nconv = sum(ritzconv)

        if nconv >= neigs
            niter = i
            break
        end

        neigsadj::Int = neigsadjusted(neigs, ncv, nconv, ritzest, prec)

        applyshifts!(neigsadj, V, H, f, ritzval)
        arnoldifac!(neigsadj, ncv, V, H, f, A, prec)
        nmatop += (ncv - neigsadj)
        ritzpairs!(which, H, ritzval, ritzvec, ritzest)
    end

    ## Final sorting of Ritz values
    ix = sortperm(ritzval[1:neigs], by = abs, rev = true)
    converged = ix[ritzconv[ix]]
    eigval = ritzval[converged]
    if returnvec
        eigvec = V * ritzvec[:, converged]
    else
        eigvec = nothing
    end

    return eigval, eigvec, nconv, niter, nmatop
end
