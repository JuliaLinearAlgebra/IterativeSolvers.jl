using Base.Test

A = [1. 2 3 4; 5 6 7 8]
S1 = svdfact(A)
S2 = rsvd(A, 2, 0)

@assert vecnorm(abs(S1[:U]) - abs(S2[:U])) <= √(eps())
@assert vecnorm(abs(S1[:Vt]) - abs(S2[:Vt])) <= √(eps())
@assert norm(S1[:S] - S2[:S]) <= √(eps())

A = [1. 2; 3 4; 5 6; 7 8]
S1 = svdfact(A)
S2 = rsvd(A, 4, 0)

@assert vecnorm(abs(S1[:U]) - abs(S2[:U])) <= √(eps())
@assert vecnorm(abs(S1[:Vt]) - abs(S2[:Vt])) <= √(eps())
@assert norm(S1[:S] - S2[:S]) <= √(eps())

A = [1. 2 3; 4 5 6; 7 8 9]
S1 = svdfact(A)
S2 = rsvd(A, 3, 0)

@assert vecnorm(abs(S1[:U]) - abs(S2[:U])) <= √(eps())
@assert vecnorm(abs(S1[:Vt]) - abs(S2[:Vt])) <= √(eps())
@assert norm(S1[:S] - S2[:S]) <= √(eps())

