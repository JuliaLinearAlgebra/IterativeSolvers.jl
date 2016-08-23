#Stationary iterative methods
#Templates, section 2.2
export jacobi, jacobi!, gauss_seidel, gauss_seidel!, sor, sor!, ssor, ssor!

jacobi(A::AbstractMatrix, b;
       tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2) =
    jacobi!(zerox(A, b), A, b; tol=tol, maxiter=maxiter)

function jacobi!(x, A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    jacobi_method!(history, x, A, b; tol=tol, maxiter=maxiter)
    x, history
end

function jacobi_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
    iter=0
    n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log)
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
        push!(log,:resnorm,resnorm)
		resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    shrink!(log)
    setmvps(log, iter)
	x
end

gauss_seidel(A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2) =
    gauss_seidel!(zerox(A, b), A, b; tol=tol, maxiter=maxiter)

function gauss_seidel!(x, A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    gauss_seidel_method!(history, x, A, b; tol=tol, maxiter=maxiter)
    x, history
end

function gauss_seidel_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
    iter=0
    n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log)
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
        push!(log,:resnorm,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    shrink!(log)
    setmvps(log, iter)
	x
end

#Successive overrelaxation
sor(A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2) =
    sor!(zerox(A, b), A, b, ω; tol=tol, maxiter=maxiter)

function sor!(x, A::AbstractMatrix, b, ω::Real;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    sor_method!(history, x, A, b, ω; tol=tol, maxiter=maxiter)
    x, history
end

function sor_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b, ω::Real;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2)
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

    iter=0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log)
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
        push!(log,:resnorm,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    shrink!(log)
    setmvps(log, iter)
	x
end

#Symmetric successive overrelaxation
#A must be symmetric
ssor(A::AbstractMatrix, b, ω::Real;
     tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)) =
    ssor!(zerox(A, b), A, b, ω; tol=tol, maxiter=maxiter)

function ssor!(x, A::AbstractMatrix, b, ω::Real;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2))
    history = ConvergenceHistory()
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    ssor_method!(history, x, A, b, ω; tol=tol, maxiter=maxiter)
    x, history
end

function ssor_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b, ω::Real;
        tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2))
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

    iter=0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log)
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
        push!(log,:resnorm,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
		copy!(xold, x)
	end
    shrink!(log)
    setmvps(log, iter)
	x
end
