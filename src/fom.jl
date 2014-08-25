#Full orthogonalization method

function fom(A, b, x₀=nothing)
	(x₀==nothing) && (x₀=convert(typeof(b), randn(length(b))))
	r₀ = b - A*x₀
	β = norm(r₀)
	v₁ = r₀ / β #Normalized form

	VH = Arnoldi(A, v₁)
	δx̂ = solve(VH)
	x₀ + β*δx̂
end

function gmres(A, b, x₀=nothing)
	(x₀==nothing) && (x₀=convert(typeof(b), randn(length(b))))
	r₀ = b - A*x₀
	β = norm(r₀)
	v₁ = r₀ / β #Normalized form

	VH = Arnoldi(A, v₁)
	δx̂ = solve₂(VH)
	x₀ + β*δx̂
end

type HessenbergMatrix{T} <: AbstractMatrix{T}
	H :: Matrix{T}
end
import Base: size

size(H::HessenbergMatrix, args...) = size(H.H, args...)

type HessenbergFact{T}
	V :: Matrix{T}
	H :: HessenbergMatrix{T}
end

HessenbergFact(V, H) = HessenbergFact(V, HessenbergMatrix(H))

function mgs(w, V)
	j=size(V,2)
	h=Array(eltype(w), j)
	for i=1:j
		v = V[:,i]
		h[i] = w ⋅ v
		w -= h[i] * v
	end
	w, h
end	

function Arnoldi{T}(A, v::Vector{T}, m::Int=length(v)) #with MGS
	V = Array(T, length(v), m)
	H = zeros(T, m+1, m)
	mend = m
	V[:,1] = v
	for j=1:m
		w, H[1:j, j] = mgs(A*V[:,j], V[:,1:j])
		nw = norm(w)
		m==j && break
		nw < length(w)*eps(T) && (mend=j; info("BREAKDOWN IN ITERATION $j"); break) #Breakdown
		V[:,j+1] = w/nw
		H[j+1,j]=nw
	end
	VH = mend==m ? HessenbergFact(V, H) : HessenbergFact(V[:,1:mend], H[1:(mend+1),1:mend])
	nVH = norm(A*VH.V - VH.V*VH.H.H[1:mend,1:mend])
	nVH < mend^2*eps(T) || warn("Error in Hessenberg factorization has norm $nVH > $(mend*eps(T))")
	VH
end

function solve(VH::HessenbergFact)
	n = size(VH.H,2)
	e₁ = zeros(eltype(VH.V), n)
	e₁[1] = 1
	y = VH.H.H[1:n,:] \ e₁ #First column of H⁻¹
	VH.V * y
end

function solve₂(VH::HessenbergFact)
	QR = qrfact(VH.H.H)
	e₁ = zeros(eltype(VH.V), size(VH.H,1))
	e₁[1] = 1
	y = QR[:R] \ (QR[:Q]'*e₁)[1:end-1]
	VH.V * y
end

function test()
	n = 140
	srand(123)
	A = float32(randn(n,n))
	b = float32(randn(n))

	for solver in [fom, gmres]
		x = solver(A, b)

		nr = norm(A*x-b)
		if nr > n*n*eps(eltype(x))
			@show A*x, b
			@show solver, nr, n*n*eps(eltype(x))
			@assert false
		end
	end
end

test()