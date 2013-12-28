export gmres, gmres!

#One Arnoldi iteration
#Optionally takes a truncation parameter l
function arnoldi!(K::KrylovSubspace; l=K.order)
    v = nextvec(K)
    n = min(length(K.v), l)
    h = zeros(eltype(v), n+1)
    v, h[1:n] = orthogonalize(v, K, n)
    h[n+1] = norm(v)
    (h[n+1]==0) && error("Arnoldi iteration terminated")
    append!(K, v/h[n+1]) #save latest Arnoldi vector, i.e. newest column of Q in the AQ=QH factorization 
    h #Return current Arnoldi coefficients, i.e. newest column of in the AQ=QH factorization
end

gmres(A, b, Pl=1, Pr=1;
      tol=sqrt(eps(typeof(real(b[1])))), maxiter::Int=1, restart::Int=min(20,length(b))) =
    gmres!(zerox(A,b), A, b, Pl, Pr; tol=tol, maxiter=maxiter, restart=restart)

function gmres!(x, A, b, Pl=1, Pr=1;
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
    H = zeros(T,n+1,restart)       #Hessenberg matrix
    w = zeros(T,n)                 #Working vector
    s = zeros(T,restart+1)         #Residual history
    J = zeros(T,restart,2)         #Givens rotation values
    a = zeros(T,restart)           #Subspace vector coefficients
    tol = tol * norm(Pl\b)         #Relative tolerance
    resnorms = zeros(typeof(real(b[1])), maxiter, restart)
    isconverged = false
    K = KrylovSubspace(x->Pl\(A*(Pr\x)), n, restart+1, T)
    for iter = 1:maxiter
        w    = Pl\(b - A*x)
        rho  = norm(w)
        s[1] = rho
        init!(K, w / rho)

        N = restart
        for j = 1:restart
            #Calculate next orthonormal basis vector in the Krylov subspace
            w = nextvec(K)
            H[1:j+1, j] = arnoldi!(K)

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
        w    = a[N] * K.v[N]

        for j = (N-1):-1:1
            a[j] = s[j]

            for k = (j+1):1:N
                a[j] -= H[j,k] * a[k]
            end

            a[j] = a[j] / H[j,j]
            w   += a[j] * K.v[j]
        end

        #Right preconditioner
        update!(x, 1, Pr\w)

        if rho < tol
            resnorms = resnorms[1:iter, :]
            isconverged = true
            break
        end
    end

    return x, ConvergenceHistory(isconverged, tol, length(resnorms), resnorms)
end

