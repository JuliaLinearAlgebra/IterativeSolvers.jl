export cg, cg!, master_cg, master_cg!

macro conjugate_gradients!(x,K,b,Pl,tol,maxiter,init, tol_check, finish)
  quote
    $init
    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = isa(Pl, Function) ? Pl(r) : Pl\r
    γ = dot(r, z)
    for iter=1:maxiter
      append!(K, p)
      q = nextvec(K)
      α = γ/dot(p, q)
      # α>=0 || throw(PosSemidefException("α=$α"))
      update!(x, α, p)
      r -= α*q
      $tol_check
      z = isa(Pl, Function) ? Pl(r) : Pl\r
      oldγ = γ
      γ = dot(r, z)
      β = γ/oldγ
      p = z + β*p
    end
    $finish
  end
end

cg(A, b, Pl=1; kwargs...) =  cg!(zerox(A,b), A, b, Pl; kwargs...)

master_cg(A, b, Pl=1; kwargs...) =  master_cg!(zerox(A,b), A, b, Pl; kwargs...)

function cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(x,K,b,Pl; tol=tol, maxiter=maxiter)
end

function master_cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    master_cg!(x,K,b,Pl; tol=tol, maxiter=maxiter)
end

function cg!(x, K::KrylovSubspace, b, Pl=1;
        tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2))
  @conjugate_gradients!(
    x,K,b,Pl,tol,maxiter,
    Void,
    (norm(r) < tol && break),
    x
  )
end

function master_cg!(x, K::KrylovSubspace, b, Pl=1;
        tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2))
  @conjugate_gradients!(
    x,K,b,Pl,tol,maxiter,
    resnorms = zeros(maxiter),  #init
    begin resnorms[iter] = norm(r)  #tol_check
      if resnorms[iter] < tol #Converged?
          resnorms = resnorms[1:iter]
          break
      end
    end,
    (x, ConvergenceHistory(0<resnorms[end]<tol, tol, K.mvps, resnorms)) #finish
  )
end
