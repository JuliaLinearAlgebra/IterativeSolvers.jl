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

function next(L::BiLanczos, S::BiLanczosState)
	⋅ = L.⋅
	v′ = L.K.A * S.v
	w′ = L.K̃.A * S.w

	α = v′⋅ S.w
	v̂ = v′ - α*S.v - S.β*S.v₋
	ŵ = w′ - α*S.w - S.δ*S.w₋
	c = v̂ ⋅ ŵ
	δ = √abs(c)
	β = c / δ
	w = ŵ / β
	v = v̂ / δ

	(α,S.β,S.δ,S.w,S.v), BiLanczosState(S.iter+1, w, S.w, v, S.v, α, β, δ)
end

done{T<:FloatingPoint}(L::BiLanczos{T}, S::BiLanczosState{T})= S.iter==length(S.w) || (S.iter>0 && abs(S.δ) < eps(T))
done{T}(L::BiLanczos{T}, S::BiLanczosState{T})= S.iter>0 && S.δ == 0

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
L = BiCG(A, b, b̃=b̃)
for it in L
	@show it
end

#Conjugate gradients
Asym = A+A'
LCG = BiCG(Asym, b, Aᵀ=Asym, b̃=A*b, innerprod=(x,y)->x⋅(A*y))
for it in LCG
	@show it
end
