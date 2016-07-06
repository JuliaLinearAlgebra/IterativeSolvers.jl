export cg, cg!, master_cg, master_cg!

cg(A, b; kwargs...) =  cg!(zerox(A,b), A, b; kwargs...)

function cg!(x, A, b; kwargs...)
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(x,K,b; kwargs...)
end

function cg!(x, K::KrylovSubspace, b; kwargs...)
    cg_method!(x,K,b; kwargs...)
    x
end

cg(::Type{Master}, A, b; kwargs...) =  cg!(Master, zerox(A,b), A, b; kwargs...)

function cg!(::Type{Master}, x, A, b; kwargs...)
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(Master, x,K,b; kwargs...)
end

function cg!(::Type{Master}, x, K::KrylovSubspace, b;
    tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2),
    verbose::Bool=false, plot=false, pl=1
    )
    log = MethodLog(maxiter)
    add!(log,:resnorm)
    cg_method!(x,K,b; pl=pl,tol=tol,maxiter=maxiter,log=log,verbose=verbose)
    shrink!(log)
    plot && showplot(log)
    x, ConvergenceHistory(isconverged(log,:resnorm,tol),tol,K.mvps,log)
end

#Make macro predicate for method functions?
function cg_method!(x,K,b;
    pl=1,tol::Real=size(K.A,2)*eps(),maxiter::Integer=size(K.A,2),
    log::MethodLog=MethodLog(),verbose::Bool=false
    )
    verbose && @printf("=== cg ===\n%4s\t%7s\n","iter","relres")
    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = pl\r
    γ = dot(r, z)
    for iter=1:maxiter
        append!(K, p)
        q = nextvec(K)
        α = γ/dot(p, q)
        # α>=0 || throw(PosSemidefException("α=$α"))
        update!(x, α, p)
        r -= α*q
        resnorm = norm(r)
        next!(log)
        push!(log,:resnorm,resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
        z = pl\r
        oldγ = γ
        γ = dot(r, z)
        β = γ/oldγ
        p = z + β*p
    end
    verbose && @printf("\n")
end
