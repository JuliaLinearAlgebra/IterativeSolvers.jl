import Base.LinAlg.BlasFloat

export eigvals_lanczos, svdvals_gkl

function lanczos!{T}(K::KrylovSubspace{T})
    m = K.n
    αs = Array(T, m)
    βs = Array(T, m-1)
    for j=1:m-1
        w = nextvec(K)
        if j>1 w -= βs[j-1]*K.v[1] end
        w, y = orthogonalize(w, K, 1)
        αs[j] = y[1]
        βs[j] = convert(T, norm(w))
        append!(K, w/βs[j])
    end
    αs[m]= dot(nextvec(K), lastvec(K))
    SymTridiagonal(αs, βs)
end

function eigvals_lanczos(A, neigs::Int=size(A,1); tol::Real=size(A,1)^3*eps(), maxiter::Int=size(A,1))
    K = KrylovSubspace(A, 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    resnorms = zeros(maxiter)
    e1 = eigvals(lanczos!(K), 1, neigs)
    for iter=1:maxiter
        e0, e1 = e1, eigvals(lanczos!(K), 1, neigs)
        resnorms[iter] = norm(e1-e0)
        if resnorms[iter] < tol
            resnorms = resnorms[1:iter]
            break
        end
    end
    e1, ConvergenceHistory(0<resnorms[end]<tol, tol, K.mvps, resnorms)
end

#Golub-Kahan-Lanczos algorithm for the singular values
#TODO duck-type for A::Function
function svdvals_gkl(A, Ap=A')
    n = size(A, 2)
    α, β = Array(eltype(A), n), Array(eltype(A), n-1)
    K = KrylovSubspace(A, 1)
    initrand!(K)
    u, v = nextvec(K), lastvec(K)
    for k=1:n
        α[k]=norm(u)
        if k==n break end
        u/=α[k]
        v=Ap*u-α[k]*v
        β[k]=norm(v)
        v/=β[k]
        u=A*v-β[k]*u
    end
    svdvals(Bidiagonal(α, β, true))
end

