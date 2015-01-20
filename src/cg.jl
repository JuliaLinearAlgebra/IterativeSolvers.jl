#######################
# Conjugate gradients #
#######################

export cg, cg!

import Base: start, next, done

@doc doc"""
Abstract iterative solver
""" ->
abstract IterativeSolver

@doc doc"""
Abstract state of iterative solver
""" ->
abstract IterationState

immutable KrylovSpace #XXX to replace KrylovSubspace
    A
    v0 :: Vector
end



@doc doc"""
Conjugate gradients iterative solver

Fields:

   `K`: Krylov space
   `t`: terminator
""" ->


immutable cg_hs <: IterativeSolver
    K :: KrylovSpace
    t :: Terminator
end

immutable cg_hs_state <: IterationState
    r :: Vector #The current residual
    p :: Vector #The current search direction
    ts:: stopcriteriastate
        #contains
        # - Squared norm of the _previous_ residual
        # - Current iteration
end

start(a::cg_hs) = cg_hs_state(a.K.v0, zeros(size(a.K.v0, 1)),
    stopcriteriastate(1, norm(a.K.v0)))

#TODO Slot in preconditioners
function next(a::cg_hs, old::cg_hs_state)
    resnorm² = dot(old.r, old.r)
    p = isfinite(old.ts.resnorm²) ? old.r + (resnorm² / old.ts.resnorm²) * old.p :
                                    old.r
    Ap = a.K.A*p
    µ = resnorm² / dot(p, Ap)
    µ*p, cg_hs_state(old.r-µ*Ap, p, stopcriteriastate(old.ts.iter+1, resnorm²))
end

done(a::cg_hs, currentstate::cg_hs_state) = done(a.t, currentstate.ts)

##################
# Generic solver #
##################

cg(A, b, Pl=1; kwargs...) = cg!(zerox(A,b), A, b, Pl; kwargs...)

# TODO return ConvergenceHistories
function cg!(x, A, b, Pl=1; tol::Real=size(A,2)*eps(), maxiter::Int=size(A,2))
    solve(A, b, x, maxiter, cg_hs)
end
