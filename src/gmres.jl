export gmres, gmres!

#One Arnoldi iteration
#Optionally takes a truncation parameter l
function arnoldi!(K::KrylovSubspace, w; l=K.order)
    v = nextvec(K)
    w = copy(v)
    n = min(length(K.v), l)
    h = zeros(eltype(v), n+1)
    v, h[1:n] = orthogonalize(v, K, n)
    h[n+1] = norm(v)
    (h[n+1]==0) && error("Arnoldi iteration terminated")
    append!(K, v/h[n+1]) #save latest Arnoldi vector, i.e. newest column of Q in the AQ=QH factorization
    h #Return current Arnoldi coefficients, i.e. newest column of in the AQ=QH factorization
end

function apply_givens!(H, J, j)
    for k = 1:(j-1)
        temp     =       J[k,1]  * H[k,j] + J[k,2] * H[k+1,j]
        H[k+1,j] = -conj(J[k,2]) * H[k,j] + J[k,1] * H[k+1,j]
        H[k,j]   = temp
    end
end

function compute_givens(a, b, i, j)
    T = eltype(a)
    p = abs(a)
    q = abs(b)

    if q == zero(T)
        return [one(T), zero(T), a]
    elseif p == zero(T)
        return [zero(T), sign(conj(b)), q]
    else
        m      = hypot(p,q)
        temp   = sign(a)
        return [p / m, temp * conj(b) / m, temp * m]
    end
end

gmres(A, b; kwargs...) = gmres!(zerox(A,b), A, b; kwargs...)

function gmres!(x, A, b; kwargs...)
    gmres_method!(x,A,b; kwargs...)
    x
end

gmres(::Type{Master}, A, b; kwargs...) = gmres!(Master, zerox(A,b), A, b; kwargs...)

function gmres!(::Type{Master}, x, A, b;
    tol=sqrt(eps(typeof(real(b[1])))),
    restart::Int=min(20,length(b)), maxiter::Int=restart,
    plot::Bool=false, kwargs...
    )
    log = ConvergenceHistory(restart)
    log[:tol] = tol
    reserve!(log,:resnorm,maxiter)
    gmres_method!(x,A,b; restart=restart,tol=tol,maxiter=maxiter,log=log, kwargs...)
    shrink!(log)
    plot && showplot(log)
    x, log
end

function gmres_method!(x, A, b;
    pl=1, pr=1, tol=sqrt(eps(typeof(real(b[1])))), restart::Int=min(20,length(b)),
    maxiter::Int=restart, verbose::Bool=false, log::MethodLog=DummyHistory()
    )
#Generalized Minimum RESidual
#Reference: http://www.netlib.org/templates/templates.pdf
#           2.3.4 Generalized Minimal Residual (GMRES)
#
#           http://www.netlib.org/lapack/lawnspdf/lawn148.pdf
#           Givens rotation based on Algorithm 1
#
#   Solve A*x=b using the Generalized Minimum RESidual Method with restarts
#
#   Effectively solves the equation inv(pl)*A*inv(pr)*y=b where x = inv(pr)*y
#
#   Required Arguments:
#       A: Linear operator
#       b: Right hand side
#
#   Named Arguments:
#       pl:       Left preconditioner
#       pr:       Right preconditioner
#       restart:  Number of iterations before restart (GMRES(restart))
#       maxiter:  Maximum total number of iterations (macroiters = ceil(maxiter/restart))
#       tol:      Convergence Tolerance
#
#   The input A (resp. pl, pr) can be a matrix, a function returning A*x,
#   (resp. inv(pl)*x, inv(pr)*x), or any type representing a linear operator
#   which implements *(A,x) (resp. \(pl,x), \(pr,x)).
    verbose && @printf("=== gmres ===\n%4s\t%4s\t%7s\n","rest","iter","relres")
    macroiter=0
    macroiters=Int(ceil(maxiter/restart))
    n = length(b)
    T = eltype(b)
    H = zeros(T,n+1,restart)       #Hessenberg matrix
    s = zeros(T,restart+1)         #Residual history
    J = zeros(T,restart,3)         #Givens rotation values
    tol = tol * norm(pl\b)         #Relative tolerance
    K = KrylovSubspace(x->pl\(A*(pr\x)), length(b), restart+1, eltype(b))
    for macroiter = 1:macroiters
        w    = pl\(b - A*x)
        s[1] = rho = norm(w)
        init!(K, w / rho)

        N = restart
        for j = 1:restart
            #Calculate next orthonormal basis vector in the Krylov subspace
            H[1:j+1, j] = arnoldi!(K, w)

            #Update QR factorization of H
            #The Q is stored as a series of Givens rotations in J
            #The R is stored in H
            #- Compute Givens rotation that zeros out bottom right entry of H
            apply_givens!(H, J, j)
            J[j,1:3] = compute_givens(H[j,j], H[j+1,j], j, j+1)
            #G = Base.LinAlg.Givens(restart+1, j, j+1, J[j,1], J[j,2], J[j,3])
            #- Zero out bottom right entry of H
            H[j,j]   = J[j,3]
            H[j+1,j] = zero(T)
            #- Apply Givens rotation j to s, given that s[j+1] = 0
            s[j+1] = -conj(J[j,2]) * s[j]
            #-conj(G.s) * s[j]
            s[j]  *= J[j,1] #G.c

            rho = abs(s[j+1])
            nextiter!(log)
            push!(log, :resnorm, rho)
            verbose && @printf("%3d\t%3d\t%1.2e\n",macroiter,j,rho)
            if (rho < tol) | (macroiter*restart+j == maxiter)
                N = j
                break
            end
        end
        verbose && @printf("\n")

        @eval a = $(VERSION < v"0.4-" ? Triangular(H[1:N, 1:N], :U) \ s[1:N] : UpperTriangular(H[1:N, 1:N]) \ s[1:N])
        w = a[1:N] * K
        update!(x, 1, pr\w) #Right preconditioner

        if (0<=rho<tol) | (macroiter*restart+N == maxiter)
            setconv(log, 0<=rho<tol)
            break
        end
    end
    setmvps(log, macroiter)
    verbose && @printf("\n")
end
