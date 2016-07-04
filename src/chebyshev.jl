export chebyshev, chebyshev!, master_chebyshev, master_chebyshev!

chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) =
    chebyshev!(zerox(A, b), A, b, λmin, λmax; kwargs...)

function chebyshev!(x, A, b, λmin::Real, λmax::Real; n::Int=size(A,2), kwargs...)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	chebyshev!(x, K, b, λmin, λmax; kwargs...)
end

function chebyshev!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real; kwargs...)
    chebyshev_method!(x, K, b, λmin, λmax; kwargs...)
    x
end

master_chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) =
    master_chebyshev!(zerox(A, b), A, b, λmin, λmax; kwargs...)

function master_chebyshev!(x, A, b, λmin::Real, λmax::Real; n::Int=size(A,2), kwargs...)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	master_chebyshev!(x, K, b, λmin, λmax; n=n, kwargs...)
end

function master_chebyshev!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
    n::Int=size(A,2), tol::Real = sqrt(eps(typeof(real(b[1])))),
    maxiter::Int = n^3, plot::Bool=false, kwargs...
    )
    log = MethodLog(maxiter)
    add!(log,:resnorm)
    chebyshev_method!(x,K,b,λmin,λmax; tol=tol,maxiter=maxiter,log=log,kwargs...)
    shrink!(log)
    plot && showplot(log)
    x, ConvergenceHistory(isconverged(log,:resnorm,tol),tol,K.mvps,log)
end

function chebyshev_method!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
    pr=1, tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = K.n^3,
    log::MethodLog=MethodLog(), verbose::Bool=false
    )
    verbose && @printf("=== chebyshev ===\n%4s\t%7s\n","iter","relres")
    local α, p
    K.order = 1
    tol = tol*norm(b)
    r = b - nextvec(K)
    d::eltype(b) = (λmax + λmin)/2
    c::eltype(b) = (λmax - λmin)/2
    for iter = 1:maxiter
    	z = pr\r
    	if iter == 1
    		p = z
    		α = 2/d
    	else
    		β = (c*α/2)^2
    		α = 1/(d - β)
    		p = z + β*p
    	end
    	append!(K, p)
    	update!(x, α, p)
    	r -= α*nextvec(K)
        resnorm = norm(r)
        next!(log)
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
	end
    verbose && @printf("\n")
end
