#Stationary iterative methods
#Templates, section 2.2
export  jacobi, jacobi!, gauss_seidel, gauss_seidel!, sor, sor!, ssor, ssor!

jacobi(A::AbstractMatrix, b; kwargs...) =
    jacobi!(zerox(A, b), A, b; kwargs...)

function jacobi!(x, A::AbstractMatrix, b; kwargs...)
    jacobi_method!(x, A, b; kwargs...)
    x
end

jacobi(::Type{Master}, A::AbstractMatrix, b; kwargs...) =
    jacobi!(Master, zerox(A, b), A, b; kwargs...)

function jacobi!(::Type{Master}, x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    jacobi_method!(x, A, b; tol=tol, log=log, maxiter=maxiter, verbose=verbose)
    shrink!(log)
    plot && showplot(log)
    x, log
end

function jacobi_method!(x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))),maxiter=size(A,2)^2,
    verbose::Bool=false, log::MethodLog=DummyHistory()
    )
    iter = 0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
		for i=1:n
			xi = z
			for j=[1:i-1;i+1:n]
				xi += A[i,j]*xold[j]
			end
			A[i,i]==0 && throw(SingularError())
			x[i]=(b[i]-xi)/A[i,i]
		end
		#check convergence
		resnorm = norm(A*x-b)
        nextiter!(log)
        push!(log,:resnorm,resnorm)
		resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    setmvps(log, iter)
end

gauss_seidel(A::AbstractMatrix, b; kwargs...) =
    gauss_seidel!(zerox(A, b), A, b; kwargs...)

function gauss_seidel!(x, A::AbstractMatrix, b; kwargs...)
    gauss_seidel_method!(x, A, b; kwargs...)
    x
end

gauss_seidel(::Type{Master}, A::AbstractMatrix, b; kwargs...) =
    gauss_seidel!(Master, zerox(A, b), A, b; kwargs...)

function gauss_seidel!(::Type{Master}, x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    gauss_seidel_method!(x, A, b; tol=tol, log=log, maxiter=maxiter, verbose=verbose)
    shrink!(log)
    plot && showplot(log)
    x, log
end

function gauss_seidel_method!(x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    verbose::Bool=false, log::MethodLog=DummyHistory()
    )
    iter = 0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
		for i=1:n
			σ=z
			for j=1:i-1
				σ+=A[i,j]*x[j]
			end
			for j=i+1:n
				σ+=A[i,j]*xold[j]
			end
			A[i,i]==0 && throw(SingularError())
			x[i]=(b[i]-σ)/A[i,i]
		end
		#check convergence
        resnorm = norm(A*x-b)
        nextiter!(log)
        push!(log,:resnorm,resnorm)
		resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    setmvps(log, iter)
end

#Successive overrelaxation
sor(A::AbstractMatrix, b, ω::Real; kwargs...) =
    sor!(zerox(A, b), A, b, ω; kwargs...)

function sor!(x, A::AbstractMatrix, b, ω::Real; kwargs...)
    sor_method!(x, A, b, ω; kwargs...)
    x
end

sor(::Type{Master}, A::AbstractMatrix, b, ω::Real; kwargs...) =
    sor!(Master, zerox(A, b), A, b, ω; kwargs...)

function sor!(::Type{Master}, x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    sor_method!(x, A, b, ω; tol=tol, log=log, maxiter=maxiter, verbose=verbose)
    shrink!(log)
    plot && showplot(log)
    x, log
end

function sor_method!(x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    verbose::Bool=false, log::MethodLog=DummyHistory()
    )
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")
    iter = 0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
		for i=1:n
			σ=z
			for j=1:i-1
				σ+=A[i,j]*x[j]
			end
			for j=i+1:n
				σ+=A[i,j]*xold[j]
			end
			A[i,i]==0 && throw(SingularError())
			σ=(b[i]-σ)/A[i,i]
			x[i]=xold[i]+ω*(σ-xold[i])
		end
		#check convergence
        resnorm = norm(A*x-b)
        nextiter!(log)
        push!(log,:resnorm,resnorm)
		resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    setmvps(log, iter)
end

#Symmetric successive overrelaxation
#A must be symmetric
ssor(A::AbstractMatrix, b, ω::Real; kwargs...) =
    ssor!(zerox(A, b), A, b, ω; kwargs...)

function ssor!(x, A::AbstractMatrix, b, ω::Real; kwargs...)
    ssor_method!(x, A, b, ω; kwargs...)
    x
end

ssor(::Type{Master}, A::AbstractMatrix, b, ω::Real; kwargs...) =
    ssor!(Master, zerox(A, b), A, b, ω; kwargs...)

function ssor!(::Type{Master}, x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2),
    plot::Bool=false, verbose::Bool=false
    )
    log = ConvergenceHistory()
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    ssor_method!(x, A, b, ω; tol=tol, log=log, maxiter=maxiter, verbose=verbose)
    shrink!(log)
    plot && showplot(log)
    x, log
end

function ssor_method!(x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2),
    verbose::Bool=false, log::MethodLog=DummyHistory()
    )
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")
    iter = 0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
		for i=1:n #Do a SOR sweep
			σ=z
			for j=1:i-1
				σ+=A[i,j]*x[j]
			end
			for j=i+1:n
				σ+=A[i,j]*xold[j]
			end
			A[i,i]==0 && throw(SingularError())
			σ=(b[i]-σ)/A[i,i]
			x[i]=xold[i]+ω*(σ-xold[i])
		end
		copy!(xold, x)
		for i=n:-1:1 #Do a backward SOR sweep
			σ=z
			for j=1:i-1
				σ+=A[i,j]*xold[j]
			end
			for j=i+1:n
				σ+=A[i,j]*x[j]
			end
			A[i,i]==0 && throw(SingularError())
			σ=(b[i]-σ)/A[i,i] #This line is missing in the Templates reference
			x[i]=xold[i]+ω*(σ-xold[i])
		end
		#check convergence
        resnorm = norm(A*x-b)
        nextiter!(log)
        push!(log,:resnorm,resnorm)
		resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    setmvps(log, iter)
end
