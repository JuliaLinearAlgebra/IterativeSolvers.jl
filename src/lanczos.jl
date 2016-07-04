import Base.LinAlg.BlasFloat

export eiglancz, master_eiglancz

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

function eiglancz(A; kwargs...)
    eigs, _ = eiglancz_method(A; kwargs...)
    eigs
end

function master_eiglancz(A;
    maxiter::Integer=size(A,1), plot::Bool=false,
    tol::Real = size(A,1)^3*eps(real(eltype(A))), kwargs...
    )
    log = MethodLog(maxiter)
    add!(log, :resnorm)
    eigs, mvps = eiglancz_method(A; maxiter=maxiter, log=log, kwargs...)
    shrink!(log)
    plot && showplot(log)
    eigs, ConvergenceHistory(isconverged(log,:resnorm,tol),tol,mvps,log)
end

function eiglancz_method(A;
    neigs::Int=size(A,1), tol::Real = size(A,1)^3*eps(real(eltype(A))),
    maxiter::Integer=size(A,1), verbose::Bool=false, log::MethodLog=MethodLog()
    )
    verbose && @printf("=== eiglancz ===\n%4s\t%7s\n","iter","relres")
    K = KrylovSubspace(A, size(A, 1), 2) #In Lanczos, only remember the last two vectors
    initrand!(K)
    e1 = eigvals(lanczos!(K), 1:neigs)
    for iter=1:maxiter
        e0, e1 = e1, eigvals(lanczos!(K), 1:neigs)
        resnorm = norm(e1-e0)
        next!(log)
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
    end
    verbose && @printf("\n")
    e1, K.mvps
end
