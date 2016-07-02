export cg, cg!, master_cg, master_cg!

cg(A, b; pl=1, pr=1, kwargs...) =  cg!(zerox(A,b), A, b; kwargs...)

function cg!(x, A, b; kwargs...)
  K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
  init!(K, x)
  cg!(x,K,b; kwargs...)
end

function cg!(x, K::KrylovSubspace, b;
  pl=1, pr=1, tol::Real=size(K.A,2)*eps(),
  maxiter::Integer=size(K.A,2), verbose::Bool=false
  )
  cg_method!(x,K,b,pl,pr;tol=tol,maxiter=maxiter,verbose=verbose)
  x
end

master_cg(A, b; kwargs...) =  master_cg!(zerox(A,b), A, b; kwargs...)

function master_cg!(x, A, b; kwargs...)
  K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
  init!(K, x)
  master_cg!(x,K,b; kwargs...)
end

function master_cg!(x, K::KrylovSubspace, b;
  pl=1, pr=1, tol::Real=size(K.A,2)*eps(),
  maxiter::Integer=size(K.A,2), verbose=false, plot=false
  )
  residuals = ResArray(maxiter)
  cg_method!(x,K,b,pl,pr; tol=tol,maxiter=maxiter,residuals=residuals,verbose=verbose)
  shrink!(residuals)
  plot && showplot(residuals)
  x, ConvergenceHistory(K.mvps,residuals,tol)
end

#Make macro predicate for method functions?
function cg_method!(x,K,b,pl,pr;
  tol::Real=size(K.A,2)*eps(),maxiter::Integer=size(K.A,2),
  residuals::Residuals=ResSingle(),verbose::Bool=true
  )
  verbose && @printf("=== cg ===\n%4s\t%7s\n","iter","relres")
  tol = tol * norm(b)
  r = b - nextvec(K)
  p = z = pl\r
  γ = dot(r, z)
  for iter=1:maxiter
    append!(K, p)
    q = nextvec(K)
    α = γ/dot(p, q)
    # α>=0 || throw(PosSemidefException("α=$α"))
    update!(x, α, p)
    r -= α*q
    push!(residuals,norm(r))
    verbose && @printf("%3d\t%1.2e\n",iter,last(residuals))
    last(residuals) < tol && break
    z = pl\r
    oldγ = γ
    γ = dot(r, z)
    β = γ/oldγ
    p = z + β*p
  end
  verbose && @printf("\n");
end
