

function qmr(A, b, x=zeros(length(b)), Pl=eye(length(b)), tol=0.001, maxit=length(b))

   r = b - A*x
   v_tld = copy(r)
    
   M1,M2 = lu(Pl)  # LU Decomposition

   y = M1 \ v_tld
   ρ = norm(y)

   w_tld = copy(r)
   z = M2' \ w_tld
   ξ = norm(z)

   γ = ε = s = d = 1   # s, d to be safe
   η = -1
   θ =  0

   normb = norm(b)   # compute outside loop
    
   for i = 1:maxit                  

      if ρ == 0 || ξ == 0
          break     # error message
      end

      v = v_tld / ρ
      y = y / ρ

      w = w_tld / ξ
      z = z / ξ

      δ = dot(z,y)
      if δ == 0
          break     #error message
      end

      y_tld = M2 \ y
      z_tld = M1'\ z

      if i == 1
         p = copy(y_tld)
         q = copy(z_tld) 
      else 
         p = y_tld - (dot(ξ,δ) / ε)*p
         q = z_tld - (dot(ρ,δ) / ε)*q
      end
       
      p_tld = A*p      # Krylov subspace, nextvec! can be used here
      ε = dot(q,p_tld)
       
      if ε == 0
          break    # error message
      end

      β = ε / δ
       
      if β == 0
          break    #error message
      end

      v_tld = p_tld - β*v
      y =  M1 \ v_tld

      ρ_1 = copy(ρ)
      ρ = norm(y)
      w_tld = (A'*q) - (β*w)
      z =  M2' \ w_tld

      ξ = norm(z)

      γ_1 = copy(γ)
      θ_1 = copy(θ)

      θ = ρ / (γ_1*β )
      γ = 1 / sqrt(1 + (θ^2))
       
      if γ == 0
          break  #error message
      end

      η = -η*ρ_1*(γ^2) / (β*(γ_1^2))

      if i == 1
         d = η*p
         s = η*p_tld 
      else
         d = η*p + (( θ_1*γ )^2)*d
         s = η*p_tld + (( θ_1*γ )^2)*s      
      end

      x = x + d                       

      r = r - s
       
      error = norm(r) / normb           
      if error <= tol
          break     # converged
      end

   end
    x
end


# random test values:

A = [1 2; 3 4]
b = [4;5]
Pl = [2 3; 4 4]
x = [1;3]

qmr(A,b)
qmr(A,b,x,Pl)
