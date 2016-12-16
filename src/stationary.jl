#Stationary iterative methods
#Templates, section 2.2
export jacobi, jacobi!, gauss_seidel, gauss_seidel!, sor, sor!, ssor, ssor!

####################
# API method calls #
####################

jacobi(A::AbstractMatrix, b; kwargs...) = jacobi!(zerox(A, b), A, b; kwargs...)

function jacobi!(x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    plot::Bool=false, log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    jacobi_method!(history, x, A, b; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

function jacobi_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))),maxiter=size(A,2)^2,
    verbose::Bool=false
    )
    verbose && @printf("=== jacobi ===\n%4s\t%7s\n","iter","resnorm")
    iter=0
    n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log,mvps=1)
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
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
		resnorm < tol && (setconv(log, resnorm>=0); break)
		@blas! xold = x
	end
    verbose && @printf("\n")
	x
end

####################
# API method calls #
####################

gauss_seidel(A::AbstractMatrix, b; kwargs...) =
    gauss_seidel!(zerox(A, b), A, b; kwargs...)

function gauss_seidel!(x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    plot::Bool=false, log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    gauss_seidel_method!(history, x, A, b; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

function gauss_seidel_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    verbose::Bool=false
    )
    verbose && @printf("=== gauss_seidel ===\n%4s\t%7s\n","iter","resnorm")
    iter=0
    n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log,mvps=1)
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
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
		@blas! xold = x
	end
    verbose && @printf("\n")
	x
end

####################
# API method calls #
####################

sor(A::AbstractMatrix, b, ω::Real; kwargs...) =
    sor!(zerox(A, b), A, b, ω; kwargs...)

function sor!(x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    plot::Bool=false, log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    sor_method!(history, x, A, b, ω; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

function sor_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2)^2,
    verbose::Bool=false
    )
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

    verbose && @printf("=== sor ===\n%4s\t%7s\n","iter","resnorm")
    iter=0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log,mvps=1)
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
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
		@blas! xold = x
	end
    verbose && @printf("\n")
	x
end

####################
# API method calls #
####################

ssor(A::AbstractMatrix, b, ω::Real; kwargs...) =
    ssor!(zerox(A, b), A, b, ω; kwargs...)

function ssor!(x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2),
    plot::Bool=false, log::Bool=false, kwargs...
    )
    (plot & !log) && error("Can't plot when log keyword is false")
    history = ConvergenceHistory(partial=!log)
    history[:tol] = tol
    reserve!(history,:resnorm, maxiter)
    ssor_method!(history, x, A, b, ω; tol=tol, maxiter=maxiter, kwargs...)
    (plot || log) && shrink!(history)
    plot && showplot(history)
    log ? (x, history) : x
end

#########################
# Method Implementation #
#########################

function ssor_method!(log::ConvergenceHistory, x, A::AbstractMatrix, b, ω::Real;
    tol=size(A,2)^3*eps(typeof(real(b[1]))), maxiter=size(A,2),
    verbose::Bool=false
    )
	0 < ω < 2 || warn("ω = $ω lies outside the range 0<ω<2 which is required for convergence")

    verbose && @printf("=== sor ===\n%4s\t%7s\n","iter","resnorm")
    iter=0
	n = size(A,2)
    xold = copy(x)
    z = zero(Amultype(A, x))
    tol = tol * norm(b)
	for iter=1:maxiter
        nextiter!(log,mvps=1)
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
		@blas! xold = x
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
        verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
        resnorm < tol && (setconv(log, resnorm>=0); break)
		@blas! xold = x
	end
    verbose && @printf("\n")
	x
end

#################
# Documentation #
#################

let
#Initialize parameters
doc1_call = """    jacobi(A, b)
"""
doc1!_call = """    jacobi!(x, A, b)
"""
doc2_call = """    gauss_seidel(A, b)
"""
doc2!_call = """    gauss_seidel!(x, A, b)
"""
doc3_call = """    sor(A, b, ω)
"""
doc3!_call = """    sor!(x, A, b, ω)
"""
doc4_call = """    ssor(A, b, ω)
"""
doc4!_call = """    ssor!(x, A, b, ω)
"""
doc1_msg = "Solve A*x=b with the Jacobi method."
doc2_msg = "Solve A*x=b with the Gauss-Seidel method."
doc3_msg = "Solve A*x=b with the successive overrelaxation method."
doc4_msg = "Solve A*x=b with the symmetric successive overrelaxation method."
doc1!_msg = "Overwrite `x`.\n\n" * doc1_msg
doc2!_msg = "Overwrite `x`.\n\n" * doc2_msg
doc3!_msg = "Overwrite `x`.\n\n" * doc3_msg
doc4!_msg = "Overwrite `x`.\n\n" * doc4_msg
doc1_arg = ""
doc2_arg = ""
doc3_arg = "`shift::Number=0`: shift to be applied to matrix A."
doc4_arg = "`shift::Number=0`: shift to be applied to matrix A."
doc1!_arg = "`x`: initial guess, overwrite final estimation."
doc2!_arg = "`x`: initial guess, overwrite final estimation."
doc3!_arg = "`x`: initial guess, overwrite final estimation.\n\n$doc3_arg"
doc4!_arg = "`x`: initial guess, overwrite final estimation.\n\n$doc4_arg"

doc1_version = (jacobi, doc1_call, doc1_msg, doc1_arg)
doc2_version = (gauss_seidel, doc2_call, doc2_msg, doc2_arg)
doc3_version = (sor, doc3_call, doc3_msg, doc3_arg)
doc4_version = (ssor, doc4_call, doc4_msg, doc4_arg)
doc1!_version = (jacobi!, doc1!_call, doc1!_msg, doc1!_arg)
doc2!_version = (gauss_seidel!, doc2!_call, doc2!_msg, doc2!_arg)
doc3!_version = (sor!, doc3!_call, doc3!_msg, doc3!_arg)
doc4!_version = (ssor!, doc4!_call, doc4!_msg, doc4!_arg)

i=0
docstring = Vector(8)

#Build docs
for (func, call, msg, arg) in [doc1_version, doc2_version, doc3_version, doc4_version,
                                doc1!_version, doc2!_version, doc3!_version, doc4!_version]
i+=1
docstring[i] = """
$call

$msg

`ch` is a `ConvergenceHistory` object. Otherwise it will only return `x`.
If `log` is set to `true` is given, method will output a tuple `x, ch`. Where
The `plot` attribute can only be used when `log` is set version.

# Arguments

$arg

`A`: linear operator.

## Keywords

`tol::Real = size(A,2)^3*eps()`: stopping tolerance.

`maxiter::Integer = size(A,2)^2`: maximum number of iterations.

`verbose::Bool = false`: verbose flag.

`log::Bool = false`: output an extra element of type `ConvergenceHistory`
containing extra information of the method execution.

`plot::Bool = false`: plot data. (Only when `log` is set)

# Output

**if `log` is `false`**

`x`: approximated solution.

**if `log` is `true`**

`x`: approximated solution.

`ch`: convergence history.

**ConvergenceHistory keys**

`:tol` => `::Real`: stopping tolerance.

`:resnom` => `::Vector`: residual norm at each iteration.

"""
end

@doc docstring[1] -> jacobi
@doc docstring[5] -> jacobi!
@doc docstring[2] -> gauss_seidel
@doc docstring[6] -> gauss_seidel!
@doc docstring[3] -> sor
@doc docstring[7] -> sor!
@doc docstring[4] -> ssor
@doc docstring[8] -> ssor!
end
