export chebyshev, chebyshev!

####################
# API method calls #
####################

"""
    chebyshev(A, b, λmin, λmax)

Solve A*x=b using the chebyshev method.

# Arguments

* `A`: linear operator.

* `b`: light hand side.

* `λmin`: minimum eigenvalue lower estimation.

* `λmax`: maximum eigenvalue upper estimation.

## Keywords

* `Pr = 1`: right preconditioner of the method.

* `tol::Real = sqrt(eps())`: stopping tolerance.

* `maxiter::Integer = size(A,2)^3`: maximum number of iterations.

* `verbose::Bool = false`: verbose flag.

# Output

* approximated solution.

"""
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

chebyshev(::Type{Master}, A, b, λmin::Real, λmax::Real; kwargs...) =
    chebyshev!(Master, zerox(A, b), A, b, λmin, λmax; kwargs...)

function chebyshev!(::Type{Master}, x, A, b, λmin::Real, λmax::Real; n::Int=size(A,2), kwargs...)
	K = KrylovSubspace(A, n, 1, Adivtype(A, b))
	init!(K, x)
	chebyshev!(Master, x, K, b, λmin, λmax; n=n, kwargs...)
end

function chebyshev!(::Type{Master}, x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
    n::Int=size(A,2), tol::Real = sqrt(eps(typeof(real(b[1])))),
    maxiter::Int = n^3, plot::Bool=false, kwargs...
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    chebyshev_method!(x,K,b,λmin,λmax; tol=tol,maxiter=maxiter,log=log,kwargs...)
    shrink!(log)
    plot && showplot(log)
    x, log
end

#########################
# Method Implementation #
#########################

function chebyshev_method!(x, K::KrylovSubspace, b, λmin::Real, λmax::Real;
    Pr=1, tol::Real = sqrt(eps(typeof(real(b[1])))), maxiter::Int = K.n^3,
    log::MethodLog=DummyHistory(), verbose::Bool=false
    )
    verbose && @printf("=== chebyshev ===\n%4s\t%7s\n","iter","relres")
    local α, p
    K.order = 1
    tol = tol*norm(b)
    r = b - nextvec(K)
    d::eltype(b) = (λmax + λmin)/2
    c::eltype(b) = (λmax - λmin)/2
    for iter = 1:maxiter
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
        resnorm = norm(r)
        nextiter!(log)
        push!(log, :resnorm, resnorm)
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && break
	end
    setmvps(log, K.mvps)
    setconv(log, 0<=norm(r)<tol)
    verbose && @printf("\n")
end
