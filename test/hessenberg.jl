using IterativeSolvers
using Base.Test

@testset "Hessenberg" begin

    # Some well-conditioned Hessenberg matrix
    H1 = [
        1.19789 1.42354 -0.0371401  0.0118481 -0.0362113  0.00269463;
        1.46142 4.01953  0.890729  -0.0157701 -0.0300656 -0.0191307;
        0.0     1.08456  3.35179    0.941966   0.0439339 -0.072888; 
        0.0     0.0      1.29071    3.1746     0.853378   0.0202058; 
        0.0     0.0      0.0        1.32227    3.06086    1.18129; 
        0.0     0.0      0.0        0.0        1.58682    2.99037;
        0.0     0.0      0.0        0.0        0.0        1.45345
    ]

    H2 = [
        1.14661+1.00222im  1.38611+0.0100307im  -0.0280623+0.00494128im    0.0205642-0.00867472im;
        1.43066+0.0im      4.11145+1.0144im       0.883563+0.00583816im  -0.00332255-0.0274278im;
            0.0+0.0im       1.0865+0.0im            3.2653+0.990342im       0.942434+0.00586481im;
            0.0+0.0im          0.0+0.0im           1.31146+0.0im             3.11491+1.01587im;
            0.0+0.0im          0.0+0.0im               0.0+0.0im             1.42175+0.0im
    ]

    for H = (H1, H2)
        T = eltype(H)

        # Fist standard basis vector as rhs
        rhs = [i == 1 ? one(T) : zero(T) for i = 1 : size(H, 1)]

        # Compare \ against the optimized version.
        solution_with_residual = copy(rhs)
        A_ldiv_B!(IterativeSolvers.FastHessenberg(copy(H)), solution_with_residual)
        solution = H \ rhs

        # First part is the solution
        @test solution_with_residual[1 : size(H, 2)] ≈ H \ rhs

        # Last element is the residual
        @test abs(last(solution_with_residual)) ≈ norm(H * solution - rhs)
    end
end