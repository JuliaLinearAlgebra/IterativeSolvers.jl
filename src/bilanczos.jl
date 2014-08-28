import Base: start, next, done

immutable K{T}
	A
	b :: AbstractVector{T}
end

#Biorthogonal Lanczos
immutable BiLanczos{T}
	K:: K{T}
	K̃:: K{T}
	⋅ :: Function
	full :: Bool #Controls output type
end

BiLanczos{T}(K::K{T}, K̃::K{T}, ⋅::Function = ⋅, full::Bool=false) =
		BiLanczos{T}(K, K̃, ⋅, full)

immutable BiLanczosState{T}
	iter :: Int
	w :: AbstractVector{T}
	w₋ :: AbstractVector{T}
	v :: AbstractVector{T}
	v₋ :: AbstractVector{T}
	α :: T
	β :: T
	δ :: T
end

immutable BiLanczosStateFull{U}
	iter :: Int
	W::AbstractMatrix{U}
	V::AbstractMatrix{U}
	T::Tridiagonal{U}
end

function start{T}(L::BiLanczos{T})
	v₁, w₁ = L.K.b, L.K̃.b
	L.full ? 
	  BiLanczosStateFull(0, reshape(w₁, length(w₁), 1),
		reshape(v₁, length(v₁), 1), Tridiagonal(T[0], T[0, 0], T[0])) : #Hack to initialize a 1x1 Tridiagonal matrix
	  BiLanczosState(0, w₁, zeros(w₁), v₁, zeros(v₁), zeros(T,3)...)
end

function next{T}(L::BiLanczos{T}, S::BiLanczosStateFull{T})
	j = S.iter
	R = j==0 ? BiLanczosState(j, S.W[:,j+1], zeros(L.K̃.b), S.V[:,j+1], zeros(L.K.b), zeros(T,3)...) :
			BiLanczosState(j, S.W[:,j+1], S.W[:,j], S.V[:,j+1], S.V[:,j], S.T[j+1,j+1], S.T[j,j+1], S.T[j+1,j])

	_, R′ = next(L, R)

	j+= 1
	if j==1
		Tri = Tridiagonal([R′.δ], [R′.α, 0], [R′.β])
	else
		S.T.d[end] = R′.α
		Tri = Tridiagonal([S.T.dl, R′.δ], [S.T.d, 0.0], [S.T.du, R′.β])
	end
	S′ = BiLanczosStateFull(j, [S.W R′.w], [S.V R′.v], Tri)
	println(j)
	println(S′.W)
	println(Tri)
	println(S′.V)
	S′, S′
end

function next(L::BiLanczos, S::BiLanczosState)
	⋅ = L.⋅
	v′ = L.K.A * S.v
	w′ = L.K̃.A * S.w

	α = v′⋅ S.w
	v̂ = v′ - S.α*S.v - S.β*S.v₋
	ŵ = w′ - S.α*S.w - S.δ*S.w₋
	c = v̂ ⋅ ŵ
	δ = √abs(c)
	β = c / δ
	w = ŵ / β
	v = v̂ / δ

	S′= BiLanczosState(S.iter+1, w, S.w, v, S.v, α, β, δ)

	S′,S′
end

done{T<:FloatingPoint}(L::BiLanczos{T}, S::BiLanczosState{T})= S.iter==length(S.w) || (S.iter>0 && abs(S.δ) < eps(T))
done{T}(L::BiLanczos{T}, S::BiLanczosState{T})= S.iter>0 && S.δ == 0

done{T<:FloatingPoint}(L::BiLanczos{T}, S::BiLanczosStateFull{T})= S.iter==size(S.W,1) || (S.iter>0 && abs(S.T[end, end-1]) < eps(T))
done{T}(L::BiLanczos{T}, S::BiLanczosStateFull{T})= S.iter>0 && S.T[end, end-1] == 0

#Biconjugate gradients

immutable BiCG{T}
	K:: K{T}
	K̃:: K{T}
	⋅ :: Function
end

BiCG(A, b; Aᵀ=A', b̃=conj(b), innerprod=⋅)=
	BiCG(K(A, b), K(Aᵀ, b̃), innerprod)

immutable BiCGState{T}
	iter :: Int
	α :: T
	β :: T
	r :: AbstractVector{T}
	r̃ :: AbstractVector{T}
	p :: AbstractVector{T}
	p̃ :: AbstractVector{T}
end

start{T}(L::BiCG{T}) = BiCGState(0, one(T), one(T), L.K.b, L.K̃.b, L.K.b, L.K̃.b)

function next(L::BiCG, S::BiCGState)
	q = L.K.A * S.p
	q̃= L.K̃.A * S.p̃
	⋅ = L.⋅

	#Biorthogonality condition
	α = S.r̃ ⋅ S.r / (S.p̃ ⋅ q)
	r = S.r - α * q
	r̃ = S.r̃ - α * q̃

	#Biconjugacy condition
	β = r̃ ⋅ r / (S.r̃ ⋅ S.r)
	p = r + β * S.p
	p̃ = r̃ + β * S.p̃

	α*p, BiCGState(S.iter+1, α, β, r, r̃, p, p̃)
end

done{T<:FloatingPoint}(L::BiCG{T}, S::BiCGState{T}) = S.iter==length(S.p) || abs(S.α)*norm(S.p) < length(S.p)^2 * eps(T)
done{T}(L::BiCG{T}, S::BiCGState{T}) = S.r == zeros(S.r)

#Tests

n=4
A=randn(n,n)
b=randn(n)
b/=norm(b)
b̃=randn(n)
b̃/=norm(b̃)

println("Raw biorthogonal Lanczos")
L = BiLanczos(K(A,b),K(A',b̃),⋅)
for it in L
	@show it
end

A = reshape([1:16],4,4)
b=[1.,0,0,0]
b̃=[1.,0,0,0]
println("Raw biorthogonal Lanczos (full)")
L = BiLanczos(K(A,b),K(A',b̃),⋅,true)
for S in L
	@show S.iter, norm(S.W[:,1:S.iter]'*L.K.A*S.V[:,1:S.iter] - full(S.T)[1:S.iter,1:S.iter])
end

println("Biconjugate gradients")
L = BiCG(A, b, b̃=b̃)
for it in L
	@show it
end

#Conjugate gradients
println("Biconjugate gradients emulating conjugate gradients on SPD A")
Asym = A'A
LCG = BiCG(Asym, b, Aᵀ=Asym, b̃=A*b, innerprod=(x,y)->x⋅(A*y))
for it in LCG
	@show it
end
