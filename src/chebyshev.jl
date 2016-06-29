export chebyshev, chebyshev!, master_chebyshev, master_chebyshev!

chebyshev(A, b, λmin::Real, λmax::Real, n::Int = size(A,2); kwargs...) =
  chebyshev!(zerox(A, b), A, b, λmin, λmax, n; kwargs...)

function chebyshev!(x, A, b, λmin::Real, λmax::Real, n::Int = size(A,2); kwargs...)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	chebyshev!(x, K, b, λmin, λmax, n; kwargs...)
end

function chebyshev!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real, n::Int = size(A,2);
	pl=1, pr=1, tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = n^3,
  verbose::Bool=false
  )
  chebyshev_method!(x, K, b, λmin, λmax;
    pl=pl,pr=pr,tol=tol,maxiter=maxiter,resnorms=zeros(1),verbose=verbose
    )
  x
end

master_chebyshev(A, b, λmin::Real, λmax::Real, n::Int = size(A,2); kwargs...) =
  chebyshev!(zerox(A, b), A, b, λmin, λmax, n; kwargs...)

function master_chebyshev!(x, A, b, λmin::Real, λmax::Real, n::Int = size(A,2); kwargs...)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	master_chebyshev!(x, K, b, λmin, λmax, n; kwargs...)
end

function master_chebyshev!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real, n::Int = size(A,2);
  pl=1, pr=1, tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = n^3,
  verbose::Bool=false, plot::Bool=false
  )
  resnorms=zeros(maxiter)
  chebyshev_method!(x,K,b,λmin,λmax;
    pl=pl,pr=pr,tol=tol,maxiter=maxiter,resnorms=resnorms,verbose=verbose)
  plot && showplot(resnorms)
  (x, ConvergenceHistory(0<resnorms[end]<tol, tol, K.mvps, resnorms)) #finish
end

function chebyshev_method!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
	pl=pl, pr=pr, tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = K.n^3,
  resnorms::Vector=zeros(1), verbose::Bool=false
  )
  verbose && @printf("=== chebyshev ===\n%4s\t%7s\n","iter","relres")
	local α, p
	K.order = 1
    tol = tol*norm(b)
	r = b - nextvec(K)
	d::eltype(b) = (λmax + λmin)/2
	c::eltype(b) = (λmax - λmin)/2
	resnorms = zeros(typeof(real(b[1])), maxiter)
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
		#Check convergence
    verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
    check(tol,resnorm,resnorms,iter) && break
	end
  verbose && @printf("\n");
	x, ConvergenceHistory(resnorms[end] < tol, tol, K.mvps, resnorms)
end
