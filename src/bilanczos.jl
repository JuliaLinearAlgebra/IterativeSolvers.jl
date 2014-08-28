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
	thin :: Bool #Controls output type
end

BiLanczos{T}(K::K{T}, K̃::K{T}, ⋅::Function = ⋅, thin::Bool=true) =
		BiLanczos{T}(K, K̃, ⋅, thin)

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

start{T}(L::BiLanczos{T}) = L.thin ? 
	BiLanczosState(0, L.K̃.b, zeros(L.K̃.b), L.K.b, zeros(L.K.b), zeros(T,3)...) :
	BiLanczosStateFull(0, reshape(L.K̃.b, length(L.K̃.b), 1),
		reshape(L.K̃.b, length(L.K̃.b), 1), Tridiagonal(T[0], T[0,0], T[0]))

function next(L::BiLanczos, S::BiLanczosStateFull)
	@show j = S.iter
	Sthin = BiLanczosState(j, S.W[:,j], S.W[:,j-1], S.V[:,j], S.V[:,j-1], T[j-1,j], T[j,j-1])
	item, thinstate = next(L, Sthin)
	Snew = BiLanczosStateFull(j+1, [S.W thinstate.w], [S.V thinstate.v],
		j==1 ? Tridiagonal(S.T.dl, S.T.d, S.T.du) : 
			   Tridiagonal([S.T.dl, Sthin.δ], [S.T.d, Sthin.α], [S.T.du, Sthin.β])
		)
	item, Snew
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

done{T<:FloatingPoint}(L::BiLanczos{T}, S::BiLanczosStateFull{T})= S.iter==length(S.w) || (S.iter>0 && abs(S.T[end, end-1]) < eps(T))
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
