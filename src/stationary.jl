#Stationary iterative methods
#Templates, section 2.2
export jacobi, gauss_seidel, sor, ssor

function jacobi(A::AbstractMatrix, b, x=nothing;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
	n = size(A,2)
	if x==nothing
		xold = zeros(eltype(b), n)
		x = Array(eltype(b), n)
	else
		xold = x
	end
	resnorms = zeros(typeof(real(b[1])), maxiter)
	for iter=1:maxiter
		for i=1:n
			x[i]=zero(eltype(b))
			for j=[1:i-1;i+1:n]
				x[i] += A[i,j]*xold[j]
			end
			A[i,i]==0 && throw(SingularError())
			x[i]=(b[i]-x[i])/A[i,i]
		end
		#check convergence
		resnorms[iter] = norm(A*x-b)
		if resnorms[iter] < tol
			resnorms=resnorms[1:iter]
			break
		end
		xold=x
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, resnorms)
end	

function gauss_seidel(A::AbstractMatrix, b, x=nothing;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
	n = size(A,2)
	if x==nothing
		xold = zeros(eltype(b), n)
		x = Array(eltype(b), n)
	else
		xold = x
	end
	resnorms = zeros(typeof(real(b[1])), maxiter)
	for iter=1:maxiter
		for i=1:n
			σ=zero(eltype(b))
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
		xold=x
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, resnorms)
end

#Successive overrelaxation
function sor(A::AbstractMatrix, b, ω::Real, x=nothing;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

	n = size(A,2)
	if x==nothing
		xold = zeros(eltype(b), n)
		x = Array(eltype(b), n)
	else
		xold = x
	end
	resnorms = zeros(typeof(real(b[1])), maxiter)
	for iter=1:maxiter
		for i=1:n
			σ=0.
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
		xold=x
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, resnorms)
end

#Symmetric successive overrelaxation
#A must be symmetric
function ssor(A::AbstractMatrix, b, ω::Real, x=nothing;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2))
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

	n = size(A,2)
	if x==nothing
		xold = zeros(eltype(b), n)
		x = Array(eltype(b), n)
	else
		xold = x
	end
	resnorms = zeros(typeof(real(b[1])), maxiter)
	for iter=1:maxiter
		for i=1:n #Do a SOR sweep
			σ=0.
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
		xold = x
		for i=n:-1:1 #Do a backward SOR sweep
			σ=0.
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
		xold = x
	end
	x, ConvergenceHistory(resnorms[end]<tol, tol, resnorms)
end