export gmres

function gmres(A, b, Pl=x->x, Pr=x->x, x=nothing;
        tol=sqrt(eps(typeof(real(b[1])))), maxiter::Int=1, restart::Int=min(20,length(b)))
#Generalized Minimum RESidual
#Reference: http://www.netlib.org/templates/templates.pdf
#           2.3.4 Generalized Minimal Residual (GMRES)
#
#           http://www.netlib.org/lapack/lawnspdf/lawn148.pdf
#           Givens rotation based on Algorithm 1
#
#   Solve A*x=b using the Generalized Minimum RESidual Method with restarts
#
#   Effectively solves the equation inv(Pl)*A*inv(Pr)*y=b where x = inv(Pr)*y
#
#   Required Arguments:
#       A: Linear operator
#       b: Right hand side
#
#   Named Arguments:
#       Pl:      Left Preconditioner
#       Pr:      Right Preconditioner
#       restart: Number of iterations before restart (GMRES(restart))
#       maxiter:  Maximum number of outer iterations
#       tol:     Convergence Tolerance
#
#   The input A (resp. Pl, Pr) can be a matrix, a function returning A*x,
#   (resp. inv(Pl)*x, inv(Pr)*x), or any type representing a linear operator
#   which implements *(A,x) (resp. \(Pl,x), \(Pr,x)).

    n = length(b)
    T = eltype(b)
    V = Array(Vector{T},restart+1) #Krylov subspace
    H = zeros(T,n+1,restart)       #Hessenberg matrix
    w = zeros(T,n)                 #Working vector
    s = zeros(T,restart+1)         #Residual history
    J = zeros(T,restart,2)         #Givens rotation values
    a = zeros(T,restart)           #Subspace vector coefficients
    if x==nothing
        x = convert(Vector{T}, randn(n))
        x /= norm(x)
    end
    A_(x) = isa(A, Function) ? A(x) : A*x 
    Pl_(x) = isa(Pl, Function) ? Pl(x) : Pl\x 
    Pr_(x) = isa(Pr, Function) ? Pr(x) : Pr\x 
    tol = tol * norm(Pl_(b))
    resnorms = zeros(typeof(real(b[1])), maxiter, restart)
    isconverged = false
    for iter = 1:maxiter
        w    = b - A_(x)
        w    = Pl_(w)
        rho  = norm(w)
        s[1] = rho
        V[1] = w / rho

        N = restart
        for j = 1:restart
            #Calculate next orthonormal basis vector in the Krylov subspace
            w = Pr_(V[j])
            w = A_(w)
            w = Pl_(w)

            #Gram-Schmidt
            for k = 1:j
                H[k,j] = dot(V[k],w)
                w     -= H[k,j] * V[k]
            end
            H[j+1,j] = norm(w)
            V[j+1]   = w / H[j+1,j]

            #Apply Givens rotation to H
            for k = 1:(j-1)
                temp     =       J[k,1]  * H[k,j] + J[k,2] * H[k+1,j]
                H[k+1,j] = -conj(J[k,2]) * H[k,j] + J[k,1] * H[k+1,j]
                H[k,j]   = temp
            end

            #Compute Givens rotation j
            p = abs(H[j,j])
            q = abs(H[j+1,j])

            if q == zero(T)
                J[j,1] = one(T)
                J[j,2] = zero(T)
            elseif p == zero(T)
                J[j,1] = zero(T)
                J[j,2] = sign(conj(H[j+1,j]))
                H[j,j] = q
            else
                m      = hypot(p,q)
                temp   = sign(H[j,j])
                J[j,1] = p / m
                J[j,2] = temp * conj(H[j+1,j]) / m
                H[j,j] = temp * m
            end

            H[j+1,j] = zero(T)

            #Apply Givens rotation j to s,(assuming s[j+1] = 0)
            s[j+1] = -conj(J[j,2]) * s[j]
            s[j]   =       J[j,1]  * s[j]

            resnorms[iter, j] = rho = abs(s[j+1])
            if rho < tol
                N = j
                break
            end
        end

        #Solve Ha=s, compute w = sum(a[i]*K.v[i])
        a[N] = s[N] / H[N,N]
        w    = a[N] * V[N]

        for j = (N-1):-1:1
            a[j] = s[j]

            for k = (j+1):1:N
                a[j] -= H[j,k] * a[k]
            end

            a[j] = a[j] / H[j,j]
            w   += a[j] * V[j]
        end

        #Right preconditioner
        x += Pr_(w)

        if rho < tol
            resnorms = resnorms[1:iter, :]
            isconverged = true
            break
        end
    end

    return x, ConvergenceHistory(isconverged, tol, resnorms) 
end

