require("../src/IterativeSolvers.jl")
using IterativeSolvers
using Base.Test

A = rand(1:10, 5, 5)
b = rand(Float32, 5)
@test IterativeSolvers.Adivtype(A, b) == Float32

#### Linear operator defined as a function
# A = cycle-back operator
function shiftback!(output, b)
    n = length(b)
    length(output) == n || error("Dimension mismatch")
    for i = 1:n-1
        output[i+1] = b[i]
    end
    output[1] = b[n]
    output
end

# A' = cycle-forward operator
function shiftfwd!(output, b)
    n = length(b)
    length(output) == n || error("Dimension mismatch")
    for i = 2:n
        output[i-1] = b[i]
    end
    output[n] = b[1]
    output
end

Atimesb = [b[end], b[1:end-1]]
Aptimesb = [b[2:end], b[1]]

A = MatrixFcn{Int}(5, 5, shiftback!)
@test eltype(A) == Int
@test size(A) == (5,5)
@test size(A,1) == 5
@test size(A,2) == 5
@test length(A) == 25
@test A*b == Atimesb
output = similar(b)
@test A_mul_B!(output, A, b) == Atimesb
@test_throws A'*b

A = MatrixCFcn{Int}(5, 5, shiftback!, shiftfwd!)
@test eltype(A) == Int
@test size(A) == (5,5)
@test size(A,1) == 5
@test size(A,2) == 5
@test length(A) == 25
@test A*b == Atimesb
output = similar(b)
@test A_mul_B!(output, A, b) == Atimesb
@test A'*b == Aptimesb
@test Ac_mul_B!(output, A, b) == Aptimesb
