function   bicgstb(A, b, tol=1e-6, maxIter=100, M1=1.0, M2=1.0)
#  x, err = bicgstb(A, b, tol, maxIter, M1, M2)
#
# bicgstab.m solves the linear system Ax=b using the
# BiConjugate Gradient Stabilized Method with preconditioning.
#
# input   A        matrix
#         x        initial guess vector
#         b        right hand side vector
#         M 1, M2  preconditioner matrices
#         maxIter   maximum number of iterations
#         tol      error tolerance
#
# output  x        solution vector
#         err      error norm

iter = 0

# initialization
 
x = zeros(length(b))
if iseltype(A,Complex); b = complex(b); end
if iseltype(b,Complex); x = complex(x); end
err = zeros(maxIter+1)
bnrm2 = norm( b )

r = b
error = norm( r ) / bnrm2; err[1] = error
alpha = 1.0
omega  = 1.0
r_tld = r
iter = 0

m1p = 0; m2p = 0;
if isa(M1,Function); m1p = 1; end
if isa(M2,Function); m2p = 1; end

ma = 0;
if isa(A,Function); ma = 1; end


for iter = 1:maxIter
    rho   = dot(r_tld,r)
    if ( rho == 0.0 )
            println("rho = 0")
            return x, err;
    end
     if ( iter > 1 )
        beta  = ( rho/rho_1 )*( alpha/omega );
        p = r + beta*( p - omega*v );
     else
        p = r;
     end

     if m1p == 0
         p_hat = M1\p
     else
         p_hat = M1(p)
     end
     if m2p == 0
         p_hat = M2\p_hat
     else
         p_hat = M2(p_hat)
     end

     v = A * p_hat;
     alpha = rho / ( dot(r_tld,v) )

     s = r - alpha*v;
     if ( norm(s)/bnrm2 < tol )
        x = x + alpha*p_hat
        resid = norm( s ) / bnrm2
        return x, err;
     end

    if m1p == 0
        s_hat = M1\s
    else
        s_hat = M1(s)
    end
    if m2p == 0
        s_hat = M2\s_hat
    else
        s_hat = M2(s_hat)
    end
    if ma == 0
        t = A*s_hat;
    elseif ma == 1
        t = A(s_hat)
    end
    omega = ( dot(t,s)) / ( dot(t,t) )
    x = x + alpha*p_hat + omega*s_hat
    r = s - omega*t
    error = norm( r ) / bnrm2
    # print(iter,"   ",error,'\n')

    err[iter+1] = error
    if ( error <= tol )
        return x, err;
    end

     if  norm(omega) < 1e-16
             return x, err
     end

     rho_1 = rho
  end

  return x, err

  end

export bicgstb
