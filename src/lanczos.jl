import Base.LinAlg.BlasFloat

export eigvals_lanczos

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

function eigvals_lanczos(A, neigs::Int=size(A,1); tol::Real = size(A,1)^3*eps(real(eltype(A))), maxiter::Integer = size(A,1))
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    e1 = eiglancz_method!(history, A, neigs; maxiter=maxiter, tol=tol)
    e1, history
end

function eiglancz_method!(log::ConvergenceHistory, A, neigs::Int=size(A,1); tol::Real = size(A,1)^3*eps(real(eltype(A))), maxiter::Integer = size(A,1))
    K = KrylovSubspace(A, size(A, 1), 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    e1 = eigvals(lanczos!(K), 1:neigs)
    for iter=1:maxiter
        nextiter!(log)
        e0, e1 = e1, eigvals(lanczos!(K), 1:neigs)
        resnorm = norm(e1-e0)
        push!(log, :resnorm, resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
    end
    shrink!(log)
    setmvps(log, K.mvps)
    e1
end
