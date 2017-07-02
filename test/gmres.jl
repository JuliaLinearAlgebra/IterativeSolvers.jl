using IterativeSolvers
using FactCheck
using Base.Test
using LinearMaps

srand(1234321)

#GMRES
facts("gmres") do

n = 10

for method = [improved_gmres, gmres]
    for T in (Float32, Float64, Complex64, Complex128)
        context("Matrix{$T}") do

        A = convert(Matrix{T}, randn(n,n))
        L = convert(Matrix{T}, randn(n,n))
        R = convert(Matrix{T}, randn(n,n))
        b = convert(Vector{T}, randn(n))
        if T <: Complex
            A += im*convert(Matrix{T}, randn(n,n))
            L += im*convert(Matrix{T}, randn(n,n))
            R += im*convert(Matrix{T}, randn(n,n))
            b += im*convert(Vector{T}, randn(n))
        end
        F = lufact(A)
        b = b/norm(b)

        # Test optimality condition: residual should be non-increasing
        x_gmres, c_gmres = method(A, b, log = true, restart = 3, maxiter = 10);
        @fact all(diff(c_gmres[:resnorm]) .<= 0.0) --> true

        # Disable this test temporarily as there is no in-place A_ldiv_B!() for standard matrices
        # x_gmres, c_gmres = method(A, b, Pl=L, Pr=R, log=true)
        # @fact c_gmres.isconverged --> true
        # @fact norm(A*x_gmres - b) --> less_than(√eps(real(one(T))))

        x_gmres, c_gmres = method(A, b, Pl=F, maxiter=1, restart=1, log=true)
        @fact c_gmres.isconverged --> true
        @fact norm(A*x_gmres - b) --> less_than(√eps(real(one(T))))

        x_gmres, c_gmres = method(A, b, Pl=Identity(), Pr=F, maxiter=1, restart=1, log=true)
        @fact c_gmres.isconverged --> true
        @fact norm(A*x_gmres - b) --> less_than(√eps(real(one(T))))
        end
    end

    for T in (Float64, Complex128)
        context("SparseMatrixCSC{$T}") do
        A = sprandn(n,n,0.5)+0.001*eye(T,n,n)
        L = sprandn(n,n,0.5)
        R = sprandn(n,n,0.5)
        b = convert(Vector{T}, randn(n))
        if T <: Complex
            A += im*sprandn(n,n,0.5)
            L += im*sprandn(n,n,0.5)
            R += im*sprandn(n,n,0.5)
            b += im*randn(n)
        end
        F = lufact(A)
        b = b / norm(b)

        # Test optimality condition: residual should be non-increasing
        x_gmres, c_gmres = method(A, b, log = true, restart = 3, maxiter = 10);
        @fact all(diff(c_gmres[:resnorm]) .<= 0.0) --> true

        # Disable this test temporarily as there is no in-place A_ldiv_B!() for standard matrices
        # x_gmres, c_gmres= method(A, b, Pl=L, Pr=R, log=true)
        # @fact c_gmres.isconverged --> true
        # @fact norm(A*x_gmres - b) --> less_than(√eps(real(one(T))))

        x_gmres, c_gmres = method(A, b, Pl=F, maxiter=1, restart=1, log=true)
        @fact c_gmres.isconverged --> true
        @fact norm(A*x_gmres - b) --> less_than(√eps(real(one(T))))

        x_gmres, c_gmres = method(A, b, Pl = Identity(), Pr=F, maxiter=1, restart=1, log=true)
        @fact c_gmres.isconverged --> true
        @fact norm(A*x_gmres - b) --> less_than(√eps(real(one(T))))
        end

        context("Linear operator defined as a function") do
            A = LinearMap(cumsum!, 100, 100, Float64; ismutating=true)
            rhs = randn(size(A,2))
            rhs/= norm(rhs)
            tol = 1e-5

            x = method(A,rhs;tol=tol,maxiter=2000)
            @fact norm(A*x - rhs) --> less_than_or_equal(tol)
        end
    end
end
end
