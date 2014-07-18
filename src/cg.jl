export cg, cg!

import Base: start, next, done
abstract IterativeSolver
abstract IterationState

immutable KrylovSpace #XXX to replace KrylovSubspace
    A
    v0 :: Vector
end

#####################
# Stopping criteria #
#####################

abstract stopcriterion

immutable stopcriteriastate{T}
    iter :: Int
    resnorm² :: T
end

immutable maxiterations <: stopcriterion
    maxiter :: Int
end
done(criterion::maxiterations, current::stopcriteriastate) = 
    current.iter>=criterion.maxiter

immutable absresnorm{T<:Real} <: stopcriterion
    threshold :: T
end
done(criterion::absresnorm, current::stopcriteriastate) =
    current.resnorm²<criterion.threshold^2

#This termination criterion gets replaced by an absresnorm
#once a resnorm is computed (its continuation is an absresnorm)
immutable relresnorm{T<:Real} <: stopcriterion
    relthreshold :: T
end
#If this criterion is still around, the starting norm is
#unknown and there's no way to know if it should stop
done(criterion::relresnorm, current::stopcriteriastate) =
    false

type Terminator
    criteria :: Vector{stopcriterion}
end
push!(T::Terminator, criterion::stopcriterion) =
    push!(T.criteria, criterion)

function done(termination::Terminator, currentstate::stopcriteriastate)
    for (idx, criterion) in enumerate(termination.criteria)
        if isa(criterion, relresnorm) && isfinite(currentstate.resnorm²)
            #There is a computed resnorm, so replace the
            #relative criterion it by its continuation as an absresnorm
            #Be careful to preserve ordering of criteria
            #since we are modifying the list while iterating over it
            deleteat!(termination.criteria, idx)
            insert!(termination.criteria, idx,
                absresnorm(criterion.relthreshold*√currentstate.resnorm²))
        else #check if current criterion is good to terminate
            done(criterion, currentstate) && return true
        end
    end
    return length(termination.criteria)==0 #If no criteria, always terminate
end

################################
# Conjugate gradients          #
# The Hestenes-Stiefel variant #
################################

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

start(a::cg_hs) = cg_hs_state(a.K.v0, zeros(size(a.K.v0, 1)), stopcriteriastate(1, norm(a.K.v0)))

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
