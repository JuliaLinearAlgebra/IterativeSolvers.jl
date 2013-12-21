export chebyshev, chebyshev!

chebyshev(A, b, λmin::Real, λmax::Real, Pr=1, n=size(A,2);
          tol::Real=sqrt(eps(typeof(real(b[1])))), maxiter::Int=n^3) =
    chebyshev!(randx(A, b), A, b, λmin, λmax, Pr, n; tol=tol,maxiter=maxiter)

function chebyshev!(x, A, b, λmin::Real, λmax::Real, Pr=1, n=size(A,2); tol::Real=sqrt(eps(typeof(real(b[1])))), maxiter::Int=n^3)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	chebyshev!(x, K, b, λmin, λmax, Pr; tol=tol, maxiter=maxiter)
end

function chebyshev!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real, Pr=1; tol::Real=sqrt(eps(typeof(real(b[1])))), maxiter::Int=K.n^3)
	K.order=1
	r = b - nextvec(K)
	d::eltype(b) = (λmax+λmin)/2
	c::eltype(b) = (λmax-λmin)/2
	resnorms=zeros(typeof(real(b[1])), maxiter)
	for iter=1:maxiter
		z = Pr\r
		if iter==1
			p = z
			α = 2/d
		else
			β = (c*α/2)^2
			α = 1/(d - β)
			p = z + β*p
		end
		append!(K, p)
		update!(x, α, p)
		r -= α * nextvec(K)
		#Check convergence
		resnorms[iter]=norm(r)
		if resnorms[iter] < tol
			resnorms = resnorms[1:iter]
			break
		end
	end
	x, ConvergenceHistory(resnorms[end] < tol, tol, K.mvps, resnorms)
end
