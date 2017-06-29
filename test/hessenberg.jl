using IterativeSolvers
using FactCheck
using LinearMaps

import IterativeSolvers: Hessenberg, solve!

facts("hessenberg") do

context("Solve") do
    # Some well-conditioned Hessenberg matrix
    H = [
        1.19789 1.42354 -0.0371401  0.0118481 -0.0362113  0.00269463;
        1.46142 4.01953  0.890729  -0.0157701 -0.0300656 -0.0191307;
        0.0     1.08456  3.35179    0.941966   0.0439339 -0.072888; 
        0.0     0.0      1.29071    3.1746     0.853378   0.0202058; 
        0.0     0.0      0.0        1.32227    3.06086    1.18129; 
        0.0     0.0      0.0        0.0        1.58682    2.99037;
        0.0     0.0      0.0        0.0        0.0        1.45345
    ]

    # Fist standard basis vector as rhs
    rhs = [i == 1 ? 1.0 : 0.0 for i = 1 : size(H, 1)]

    # Compare \ against the optimized version.
    solution_with_residual = copy(rhs)
    IterativeSolvers.solve!(IterativeSolvers.Hessenberg(copy(H)), solution_with_residual)
    solution = H \ rhs

    # First part is the solution
    @fact solution_with_residual[1 : size(H, 2)] --> roughly(H \ rhs)

    # Last element is the residual
    @fact last(solution_with_residual) --> roughly(norm(H * solution - rhs))
end

end