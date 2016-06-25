using UnicodePlots
import Base: \
export cg, cg!, master_cg, master_cg!, conjugate_gradients!

\(f::Function, b::Vector) = f(b)

#@printf("%3d\t%1.2e\n",iter,resvec[iter])

#Borders are always white, issue on UnicodePlots? light terminals suffer
function showplot(vals::Vector)
  isdefined(Main, :UnicodePlots) || warn("UnicodePlots not found; no plotsies T.T ")
  println(lineplot(1:length(vals), vals, title = "Convergence", name = "resnorm"))
end

check(tol::Real, resnorm::Real, ::Vector, ::Integer, ::Type{Val{true}}) = resnorm < tol

function check(tol::Real, resnorm::Real, resnorms::Vector, iter::Int, ::Type{Val{false}})
  resnorms[iter] = resnorm
  if resnorms[iter] < tol
    resize!(resnorms,iter)
    return true
  end
  false
end

check(tol::Real, resnorm::Real, resnorms::Vector, iter::Integer) =
  check(tol,resnorm,resnorms,iter,Val{length(resnorms)==1})

#Make macro predicate for method functions?
function conjugate_gradients!(x,K,b,Pl,tol,maxiter,resnorms,verbose)
  println("=== cg ===")
  @printf("%4s\t%7s\n","iter","relres")
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
    resnorm = norm(r)
    verbose && @printf("%3d\t%1.2e\n",iter,resnorm)
    check(tol,resnorm,resnorms,iter) && break
    z = Pl\r
    oldγ = γ
    γ = dot(r, z)
    β = γ/oldγ
    p = z + β*p
  end
  verbose && @printf("\n");
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
