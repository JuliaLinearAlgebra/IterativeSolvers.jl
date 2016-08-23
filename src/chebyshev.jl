export chebyshev, chebyshev!

chebyshev(A, b, λmin::Real, λmax::Real, Pr = 1, n = size(A,2);
          tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = n^3) =
    chebyshev!(zerox(A, b), A, b, λmin, λmax, Pr, n; tol=tol,maxiter=maxiter)

function chebyshev!(x, A, b, λmin::Real, λmax::Real, Pr = 1, n = size(A,2);
	tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = n^3
    )
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm,maxiter)

	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	chebyshev_method!(history, x, K, b, λmin, λmax, Pr; tol = tol, maxiter = maxiter)
    x, history
end

function chebyshev_method!(
    log::ConvergenceHistory, x, K::KrylovSubspace, b, λmin::Real, λmax::Real,
    Pr = 1; tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = K.n^3
    )

	local α, p
	K.order = 1
    tol = tol*norm(b)
	r = b - nextvec(K)
	d::eltype(b) = (λmax + λmin)/2
	c::eltype(b) = (λmax - λmin)/2
	for iter = 1:maxiter
        nextiter!(log)
		z = Pr\r
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
		#Check convergence
		resnorm = norm(r)
        push!(log, :resnorm, resnorm)
        resnorm < tol && break
	end
    shrink!(log)
    setmvps(log, K.mvps)
    setconv(log, 0<=norm(r)<tol)
	x
end
