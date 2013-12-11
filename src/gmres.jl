export gmres

function gmres{T<:BlasFloat}(M1::Function, A::Function, M2::Function, b::Vector{T}; x::Array{T,1} = zeros(T,length(b)), restart::Int = min(20,length(b)), max_it::Int = 1, tol = 1e-6)
#Generalized Minimum RESidual
#Reference: http://www.netlib.org/templates/templates.pdf
#           2.3.4 Generalized Minimal Residual (GMRES)
#
#           http://www.netlib.org/lapack/lawnspdf/lawn148.pdf
#           Givens rotation based on Algorithm 1
#
#   Solve A*x=b using the Generalized Minimum RESidual Method with restarts
#
#   Effectively solves the equation inv(M1)*A*inv(M2)*y=b where x = inv(M2)*y
#
#   Required Arguments:
#       A: Linear operator
#       b: Right hand side
#
#   Named Arguments:
#       M1:      Left Preconditioner
#       M2:      Right Preconditioner
#       restart: Number of iterations before restart (GMRES(restart))
#       max_it:  Maximum number of outer iterations
#       tol:     Convergence Tolerance
#
#   The input A (resp. M1, M2) can be a matrix, a function returning A*x,
#   (resp. inv(M1)*x, inv(M2)*x), or any type representing a linear operator
#   which implements *(A,x) (resp. \(M1,x), \(M2,x)).

    n = length(b)
    V = Array(Vector{T},restart+1) #Krylov subspace
    H = zeros(T,n+1,restart)       #Hessenberg matrix
    w = zeros(T,n)                 #Working vector
    s = zeros(T,restart+1)         #Residual history
    J = zeros(T,restart,2)         #Givens rotation values
    a = zeros(T,restart)           #Subspace vector coefficients

    tol = tol * norm(M1(b))
    for i = 1:max_it
        w    = b - A(x)
        w    = M1(w)
        rho  = norm(w)
        s[1] = rho
        V[1] = w / rho

        N = restart
        for j = 1:restart
            #Calculate next orthonormal basis vector in the Krylov subspace
            w = M2(V[j])
            w = A(w)
            w = M1(w)

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

            rho = abs(s[j+1])
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
        x += M2(w)

        if rho < tol
            break
        end
    end

    return x
end
gmres(A, b;M1 = (x->x), M2 = (x->x), args...) = gmres(M1,       A,       M2,       b;args...)
gmres(M1,          A,          M2, b;args...) = gmres(x->(M1\x),A,       M2,       b;args...)
gmres(M1::Function,A,          M2, b;args...) = gmres(M1,       x->(A*x),M2,       b;args...)
gmres(M1::Function,A::Function,M2, b;args...) = gmres(M1,       A,       x->(M2\x),b;args...)
