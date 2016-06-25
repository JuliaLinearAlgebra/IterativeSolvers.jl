using UnicodePlots
import Base: \
export cg, cg!, master_cg, master_cg!, conjugate_gradients!

\(f::Function, b::Vector) = f(b)

#Borders are always white, issue on ubnicode plots? light terminals suffer
function showplot(vals::Vector)
  isdefined(Main, :UnicodePlots) || warn("UnicodePlots not found; no plotsies T.T ")
  println(lineplot(1:length(vals), vals, title = "Convergence", name = "resnorm"))
end

check(tol::Real, r::Vector, ::Vector, ::Integer, ::Type{Val{true}}) = norm(r) < tol

function check(tol::Real, r::Vector, resnorms::Vector, iter::Int, ::Type{Val{false}})
  resnorms[iter] = norm(r)
  if resnorms[iter] < tol
      resize!(resnorms,iter)
      return true
  end
  false
end

check(tol::Real, r::Vector, resnorms::Vector, iter::Integer) =
  check(tol,r,resnorms,iter,Val{length(resnorms)==1})

#Make macro predicate for method functions?
function conjugate_gradients!(x,K,b,Pl,tol,maxiter,resnorms,verbose)
    tol = tol * norm(b)
    r = b - nextvec(K)
    p = z = Pl\r
    γ = dot(r, z)
    for iter=1:maxiter
      append!(K, p)
      q = nextvec(K)
      α = γ/dot(p, q)
      # α>=0 || throw(PosSemidefException("α=$α"))
      LinAlg.axpy!(α, p, x)   #Is it ok in performance for this function to automatically promote?
      r -= α*q
      #showStatus()
      check(tol,r,resnorms,iter) && break
      z = Pl\r
      oldγ = γ
      γ = dot(r, z)
      β = γ/oldγ
      p = z + β*p
    end
end

cg(A, b, Pl=1; kwargs...) =  cg!(zerox(A,b), A, b, Pl; kwargs...)

master_cg(A, b, Pl=1; kwargs...) =  master_cg!(zerox(A,b), A, b, Pl; kwargs...)

function cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2), verbose=false, plot=false)
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    cg!(x,K,b,Pl; tol=tol, maxiter=maxiter, verbose=verbose)
end

function master_cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2), verbose=false, plot=false)
    K = KrylovSubspace(A, length(b), 1, Vector{Adivtype(A,b)}[])
    init!(K, x)
    master_cg!(x,K,b,Pl; tol=tol, maxiter=maxiter, verbose=verbose, plot=plot)
end

function cg!(x, K::KrylovSubspace, b, Pl=1;
        tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2), verbose=false)
  conjugate_gradients!(x,K,b,Pl,tol,maxiter,zeros(1),verbose)
  x
end

function master_cg!(x, K::KrylovSubspace, b, Pl=1;
        tol::Real=size(K.A,2)*eps(), maxiter::Integer=size(K.A,2), verbose=false, plot=false)
  resnorms=zeros(maxiter)
  conjugate_gradients!(x,K,b,Pl,tol,maxiter,resnorms,verbose)
  plot && showplot(resnorms)
  (x, ConvergenceHistory(0<resnorms[end]<tol, tol, K.mvps, resnorms)) #finish
end
