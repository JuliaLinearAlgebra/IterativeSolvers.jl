module LinearSystemsBench

import Base.A_ldiv_B!, Base.\

using BenchmarkTools
using IterativeSolvers

include("../test/advection_diffusion.jl")

# A DiagonalMatrix that doesn't check whether it is singular in the \ op.
immutable DiagonalPreconditioner{T}
    diag::Vector{T}
end

function A_ldiv_B!{T}(y::AbstractVector{T}, A::DiagonalPreconditioner{T}, b::AbstractVector{T})
    for i = 1 : length(b)
        @inbounds y[i] = A.diag[i] \ b[i]
    end
    y
end

(\)(D::DiagonalPreconditioner, b::AbstractVector) = D.diag .\ b

function posdef(n)
    A = SymTridiagonal(fill(2.01, n), fill(-1.0, n))
    b = A * ones(n)
    return A, b
end

function cg(; n = 1_000_000, tol = 1e-6, maxiter::Int = 200)
    A, b = posdef(n)
    P = DiagonalPreconditioner(collect(linspace(1.0, 2.0, n)))

    println("Symmetric positive definite matrix of size ", n)
    println("Eigenvalues in interval [0.01, 4.01]")
    println("Tolerance = ", tol, "; max #iterations = ", maxiter)
    
    # Dry run
    initial = rand(n)
    IterativeSolvers.cg!(copy(initial), A, b, Pl = P, maxiter = maxiter, tol = tol, log = false)

    # Actual benchmark
    @benchmark IterativeSolvers.cg!(x0, $A, $b, Pl = $P, maxiter = $maxiter, tol = $tol, log = false) setup=(x0 = copy($initial))
end

function indefinite(n)
    # Generate an indefinite "hard" matrix
    srand(1)
    A = speye(n) + sprand(n, n, 5.0 / n)
    A = (A + A') / 2
    x = ones(n)
    b = A * x

    A, b
end

function gmres(; n = 100_000, tol = 1e-5, restart::Int = 15, maxiter::Int = 210)
    A, b = indefinite(n)

    println("Matrix of size ", n, " with ~", nnz(A) / n, " nonzeros per row")
    println("Tolerance = ", tol, "; restart = ", restart, "; max #iterations = ", maxiter)

    IterativeSolvers.gmres(A, b, tol = tol, restart = restart, maxiter = maxiter, verbose = true)

    @benchmark IterativeSolvers.gmres($A, $b, tol = $tol, restart = $restart, maxiter = $maxiter)
end

function bicgstabl()
    A, b = advection_dominated()

    b1 = @benchmark IterativeSolvers.bicgstabl($A, $b, 2, max_mv_products = 1000, convex_combination = false) setup = (srand(1))
    b2 = @benchmark IterativeSolvers.bicgstabl($A, $b, 2, max_mv_products = 1000, convex_combination = true) setup = (srand(1))
    b3 = @benchmark IterativeSolvers.bicgstabl($A, $b, 4, max_mv_products = 1000, convex_combination = false) setup = (srand(1))
    b4 = @benchmark IterativeSolvers.bicgstabl($A, $b, 4, max_mv_products = 1000, convex_combination = true) setup = (srand(1))

    b1, b2, b3, b4
end

function minres(n = 100_000)
    A = SymTridiagonal(fill(2.1, n), fill(-1.0, n))
    x = ones(n)
    b = A * x

    @benchmark IterativeSolvers.minres($A, $b, maxiter = 100)
end

end
