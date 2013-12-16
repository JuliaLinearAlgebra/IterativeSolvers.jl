function qmr(A,b,n,maxiter,Pl=M1,Pr=M2) #without lookahead
	#Step 0. Initialize
	xold=zeros(n) #guess
	r = b - A*xold
	vt=r; y = Pl\vt; ρ[1]=norm(y)
	wt=r #arbitrary
	z = Pr'\wt; ξ[1] = norm(z)
	γ0 = 1.; η0 = -1.
	for i=1:maxiter
		#Step 1. Generate Lanczos vectors v[], w[]
		(ρ[i]==0 || ξ[i]==0) && error("fail")
		v[i]=vt[i]/ρ[i]; y/=ρ[i]
		w[i]=wt[i]/ξ[i]; z/=ξ[i]
		δ[i]=dot(z,y); δ[i]==0 && error("fail")
		yt = Pr\y; zt = Pl'\z
		if i==1
			p[1]==yt
			q[1]==zt
		else
			p[i]=yt - (ξ[i]*δ[i]/ε[i−1]) * p[i-1]
			q[i]=yt - (ρ[i]*δ[i]/ε[i−1]) * q[i-1]
		end
		pt = A*p[i]
		ξ[i] = dot(q[i], pt); ξ[i]==0 && error("fail")
		β[i] = ξ[i]/δ[i]; β[i]==0 && error("fail")
		vt[i+1] = pt - β[i]*v[i]
		y = Pl \ vt[i+1]
		ρ[i+1]=norm(y)
		wt[i+1]=A'*q[i]-β[i]*w[i]
		z = Pr' \ wt[i+1]
		ξ[i+1] = norm(z)
		#Step 2. Update QR factorization
		θ[i] = ρ[i]+1/(γ[i]−1*abs(β[i])); γ[i] = 1/sqrt(1. + θ[i]^2); γ[i]==0 && error("fail")
		η[i] = −η[i−1]*ρ[i]*γ[i]^2/(β[i]*γ[i-1]^2)
		if i == 1
			d[1] = η[1]*p[1]; s[1] = η[1]*pt
		else
			d[i] = η[i]*p[i] + (θ[i−1]*γ[i])^2*d[i−1]
			s[i] = η[i]*pt+ (θ[i-1]*γ[i])^2*s[i−1]
		end
		#Step 3. Update solution and residual
		x[i] = x[i−1] + d[i]
		r[i] = r[i−1] − s[i]
		#Step 4. Check convergence; continue if necessary
	end
end
