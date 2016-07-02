#Stationary iterative methods
#Templates, section 2.2
export  jacobi, jacobi!, master_jacobi, master_jacobi!,
        gauss_seidel, gauss_seidel!, master_gauss_seidel, master_gauss_seidel!
        sor, sor!, master_sor, master_sor!
        ssor, ssor!, master_ssor, master_ssor!

jacobi(A::AbstractMatrix, b; kwargs...) =
    jacobi!(zerox(A, b), A, b; kwargs...)

function jacobi!(x,A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    verbose::Bool=false
    )
    jacobi_method!(x, A, b; tol=tol, maxiter=maxiter, verbose=verbose)
    x
end

master_jacobi(A::AbstractMatrix, b; kwargs...) =
    jacobi!(zerox(A, b), A, b; kwargs...)

function master_jacobi!(x,A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    verbose::Bool=false, plot::Bool=false
    )
    rest_resarray=RestResArray(maxiter, restart)
    K = KrylovSubspace(x->pl\(A*(pr\x)), length(b), restart+1, eltype(b))
    jacobi_method!(x, A, b;
        tol=tol, residuals=ResArray(maxiter), maxiter=maxiter, verbose=verbose
        )
    resnorms = extract!(rest_resarray)
    plot && showplot(resnorms)
    x, ConvergenceHistory(0<=resnorms[end]<tol, tol, length(resnorms), resnorms)
end

function jacobi_method!(x, A::AbstractMatrix, b;
    residuals::Residuals=ResSingle(), tol=size(A,2)^3*eps(typeof(real(b[1]))),
    maxiter=size(A,2)^2, verbose::Bool=false
    )
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	resnorms = zeros(typeof(real(b[1])), maxiter)
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
		resnorms[iter] = norm(A*x-b)
		if resnorms[iter] < tol
			resnorms=resnorms[1:iter]
			break
		end
		copy!(xold, x)
	end
end

gauss_seidel(A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2) =
    gauss_seidel!(zerox(A, b), A, b; tol=tol, maxiter=maxiter)

function gauss_seidel!(x, A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	resnorms = zeros(typeof(real(b[1])), maxiter)
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
		resnorms[iter] = norm(A*x-b)
		if resnorms[iter] < tol
			resnorms=resnorms[1:iter]
			break
		end
		copy!(xold, x)
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, length(resnorms), resnorms)
end

#Successive overrelaxation
sor(A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2) =
    sor!(zerox(A, b), A, b, ω; tol=tol, maxiter=maxiter)

function sor!(x, A::AbstractMatrix, b, ω::Real;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	resnorms = zeros(typeof(real(b[1])), maxiter)
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
		resnorms[iter] = norm(A*x-b)
		if resnorms[iter] < tol
			resnorms=resnorms[1:iter]
			break
		end
		copy!(xold, x)
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, length(resnorms), resnorms)
end

#Symmetric successive overrelaxation
#A must be symmetric
ssor(A::AbstractMatrix, b, ω::Real;
     tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)) =
    ssor!(zerox(A, b), A, b, ω; tol=tol, maxiter=maxiter)

function ssor!(x, A::AbstractMatrix, b, ω::Real;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2))
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	resnorms = zeros(typeof(real(b[1])), maxiter)
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
		resnorms[iter] = norm(A*x-b)
		if resnorms[iter] < tol
			resnorms=resnorms[1:iter]
			break
		end
		copy!(xold, x)
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, length(resnorms), resnorms)
end
