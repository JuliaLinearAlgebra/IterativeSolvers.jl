immutable HessenbergQR{T}
    cosθ::Vector{T}
    sinθ::Vector{T}
end

## QR decomposition on an upper Hessenberg matrix
## H will be overwritten by Q'HQ=RQ
function hessenqr!{T}(H::Matrix{T})
    ## Size of matrix
    n = size(H, 1)
    if n != size(H, 2)
        error("matrix must be square")
    end

    ## Rotation factors
    ## Each pair (cosθ[i], sinθ[i]) forms a rotation matrix
    ## Gi = [cosθ[i]  -sinθ[i]]
    ##      [sinθ[i]   cosθ[i]]
    cosθ = ones(T, n)
    sinθ = zeros(T, n)

    for i = 1:(n - 1)
        ## Make sure H is upper Hessenberg
        ## Zero the elements below H[i+1, i]
        H[(i + 2):end, i] = zero(T)
        ## Calculate cosθ and sinθ
        x = H[i, i]
        y = H[i + 1, i]
        r = hypot(x, y)
        ## If r is too small, (cosθ, sinθ) stores the original values (1, 0)
        if r < eps(T)
            r = zero(T)
        else
            cosθ[i] = x / r
            sinθ[i] = -y / r
        end
        ## Apply the rotation on the left H -> Gi * H
        ## H[i, :]     <- cosθ[i] * H[i, :] - sinθ[i] * H[i +1, :]
        ## H[i + 1, :] <- sinθ[i] * H[i, :] + cosθ[i] * H[i +1, :]
        H[i, i] = r
        H[i + 1, i] = zero(T)
        c = cosθ[i]
        s = sinθ[i]
        for j = (i + 1):n
            tmp = H[i, j]
            H[i, j]     = c * tmp - s * H[i + 1, j]
            H[i + 1, j] = s * tmp + c * H[i + 1, j]
        end
    end

    ## Apply the rotations on the right H -> H * Gi'
    ## H[:, i]     <- cosθ[i] * H[:, i] - sinθ[i] * H[:, i + 1]
    ## H[:, i + 1] <- sinθ[i] * H[:, i] + cosθ[i] * H[:, i + 1]
    for i = 1:(n - 1)
        for j = 1:(i + 1)
            tmp = H[j, i]
            H[j, i]     = cosθ[i] * tmp - sinθ[i] * H[j, i + 1]
            H[j, i + 1] = sinθ[i] * tmp + cosθ[i] * H[j, i + 1]
        end
    end

    return HessenbergQR(cosθ, sinθ)
end

## QR decomposition on a tridiagonal matrix
## H will be overwritten by Q'HQ=RQ
function tridiagqr!{T}(H::Matrix{T})
    ## Size of matrix
    n = size(H, 1)
    if n != size(H, 2)
        error("matrix must be square")
    end

    ## Force H to be tridiagonal
    maindiag = diag(H)
    subdiag = diag(H, -1)
    fill!(H, zero(T))
    for i = 1:(n - 1)
        H[i, i] = maindiag[i]
        H[i + 1, i] = subdiag[i]
        H[i, i + 1] = subdiag[i]
    end
    H[n, n] = maindiag[n]

    ## Rotation factors
    ## Each pair (cosθ[i], sinθ[i]) forms a rotation matrix
    ## Gi = [cosθ[i]  -sinθ[i]]
    ##      [sinθ[i]   cosθ[i]]
    cosθ = ones(T, n)
    sinθ = zeros(T, n)

    for i = 1:(n - 1)
        ## Calculate cosθ and sinθ
        x = H[i, i]
        y = H[i + 1, i]
        r = hypot(x, y)
        ## If r is too small, (cosθ, sinθ) stores the original values (1, 0)
        if r < eps(T)
            r = zero(T)
        else
            cosθ[i] = x / r
            sinθ[i] = -y / r
        end
        ## Apply the rotation on the left H -> Gi * H
        ## Update H[i, i] and H[i + 1, i]
        H[i, i] = r
        H[i + 1, i] = zero(T)
        ## Update H[i, i + 1] and H[i + 1, i + 1]
        c = cosθ[i]
        s = sinθ[i]
        tmp = H[i, i + 1]
        H[i,     i + 1] = c * tmp - s * H[i + 1, i + 1]
        H[i + 1, i + 1] = s * tmp + c * H[i + 1, i + 1]
        ## Update H[i, i + 2] and H[i + 1, i + 2]
        if i < n - 1
            H[i,     i + 2] = -s * H[i + 1, i + 2]
            H[i + 1, i + 2] *= c
        end
    end

    ## Apply the rotations on the right H -> H * Gi'
    ## H[:, i]     <- cosθ[i] * H[:, i] - sinθ[i] * H[:, i + 1]
    ## H[:, i + 1] <- sinθ[i] * H[:, i] + cosθ[i] * H[:, i + 1]
    for i = 1:(n - 1)
        c = cosθ[i]
        s = sinθ[i]

        tmp = H[i + 1, i]
        H[i, i]         = c * H[i, i] - s * H[i, i + 1]
        H[i + 1, i]     = c * tmp     - s * H[i + 1, i + 1]
        H[i + 1, i + 1] = s * tmp     + c * H[i + 1, i + 1]

        H[i, i + 1] = H[i + 1, i]
        if i < n - 1
            H[i, i + 2] = zero(T)
        end
    end

    return HessenbergQR(cosθ, sinθ)
end

## Apply the QR factorization to the right of a matrix A
## A -> A * Q'
function applyright!{T}(qr::HessenbergQR{T}, A::Matrix{T})
    cosθ = qr.cosθ
    sinθ = qr.sinθ
    n = length(cosθ)
    ## A[:, i]     <- cosθ[i] * A[:, i] - sinθ[i] * A[:, i + 1]
    ## A[:, i + 1] <- sinθ[i] * A[:, i] + cosθ[i] * A[:, i + 1]
    for i = 1:(n - 1)
        for j = 1:size(A, 1)
            tmp = A[j, i]
            A[j, i]     = cosθ[i] * tmp - sinθ[i] * A[j, i + 1]
            A[j, i + 1] = sinθ[i] * tmp + cosθ[i] * A[j, i + 1]
        end
    end
end
