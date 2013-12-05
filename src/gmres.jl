function gmres( A, b, tol = 1e-6, max_it = 100, restrt = 10, M = 1.0)
#
#
#
# (x, eror, iter, flag) = gmres( A, b, restrt, max_it, tol, M, x )
#
# gmres solves the linear system Ax=b
# using the Generalized Minimal residual ( GMRESm ) method with restarts .
#
# input   A        REAL nonsymmetric positive definite matrix
#         x        REAL initial guess vector
#         b        REAL right hand side vector
#         M        REAL preconditioner matrix
#         restrt   INTEGER number of iterations between restarts
#         max_it   INTEGER maximum number of iterations
#         tol      REAL eror tolerance
#
# output  x        REAL solution vector
#         eror    REAL eror norm
#         iter     INTEGER number of iterations performed
#         flag     INTEGER: 0 = solution found to tolerance
#                           1 = no convergence given max_it

# initialization
iter = 0
flag = 0

x = zeros(size(b))
fm = 0; if isa(M,Function); fm = 1; end
fa = 0; if isa(A,Function); fa = 1; end
if isa(A[1,1],Complex); b = complex(b); x = complex(x); fc = 1; end

bnrm2 = norm(b)
if bnrm2 == 0.0; bnrm2 = 1.0; end

if fm == 1
    r = M(b)
else
    r = M\b
end

eror = norm( r ) / bnrm2
if eror < tol;
   return x, eror
end

n  = length(b)
V  = zeros(n,restrt+1); if fc == 1; V = complex(V); end
H  = zeros(restrt+1,restrt);  if fc == 1; H = complex(H); end
cs = zeros(restrt);  if fc == 1; cs = complex(cs); end
sn = zeros(restrt);   if fc == 1; sn = complex(sn); end
e1 = zeros(n)
e1[1] = 1.0
err = zeros(restrt*max_it)

cnt = 1
# begin iteration
for iter = 1:max_it
    if fa == 1
        r = b - A(x)
    else
        r = b - A*x
    end
    if fm == 1
        r = M(r)
    else
        r = M\r
    end

    V[:,1] = r / norm( r )
    s = norm( r )*e1;  if fc == 1; s = complex(s); end
    println("------ iter ", iter, "---------")
    for i = 1:restrt
        if fa == 1
            w = A(V[:,i])
        else
            w = A*V[:,i]
        end
        if fm == 1
            w = M(w);
        else
            w = M\w;
        end
        # basis using Gram-Schmidt
        for k = 1:i
           H[k,i] = dot(w,V[:,k])
           w = w - H[k,i]*V[:,k]
        end
        H[i+1,i] = norm( w )
        V[:,i+1] = w / H[i+1,i]
        # apply Givens rotation
        for k = 1:i-1
           temp     =  cs[k]*H[k,i] + sn[k]*H[k+1,i]
           H[k+1,i] = -sn[k]*H[k,i] + cs[k]*H[k+1,i]
           H[k,i]   = temp
        end

        # Approximate residual norm
        (cs[i],sn[i]) = rotmat( H[i,i], H[i+1,i] )
        temp   = cs[i]*s[i]
        s[i+1] = -sn[i]*s[i]
        s[i]   = temp
        H[i,i] = cs[i]*H[i,i] + sn[i]*H[i+1,i]
        H[i+1,i] = 0.0
        eror  = abs(s[i+1]) / bnrm2
        println(iter,"    ",i,"     ", eror)
        err[cnt] = eror
        cnt = cnt+1

        if eror <= tol
            y = H[1:i,1:i] \ s[1:i]
            x = x + V[:,1:i]*y
            flag = 1
            return x, err, flag
        end
    end
    if  eror <= tol;
        flag = 1
        return x, err, flag
    end
    y = H[1:restrt,1:restrt] \ s[1:restrt]
    x = x + V[:,1:restrt]*y
    if fa == 1
        r = b - A(x)
    else
        r = b - A*x
    end
    if fm == 1
        r = M(r)
    else
        r = M\r
    end
    s[restrt+1] = norm(r)
    eror = abs(s[restrt+1]) / bnrm2;

    if eror <= tol;
       flag = 1
       return x, err, flag
    end;
end

if eror > tol; flag = 1; end;

return x, err, flag

end
# end gmres.m


function rotmat( a, b )
#  [c,s] = rotmat(a,b)
# Givens rotation matrix parameters for a and b.

if b == 0.0
 c = 1.0
 s = 0.0
elseif abs(b) > abs(a)
 temp = a / b
 s = 1.0 / sqrt( 1.0 + temp^2 )
 c = temp * s
else
 temp = b / a
 c = 1.0 / sqrt( 1.0 + temp^2 )
 s = temp * c
end

return c, s

end

export gmres
