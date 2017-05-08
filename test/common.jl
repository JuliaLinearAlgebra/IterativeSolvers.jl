using IterativeSolvers
using FactCheck
using LinearMaps

facts("Basic operations") do

context("Adivtype") do
    A = rand(1:10, 5, 5)
    b = rand(Float32, 5)
    @fact IterativeSolvers.Adivtype(A, b) --> Float32
end

context("Linear operator defined as a function") do
    # A = cycle-back operator
    function shiftback!(output, b)
        n = length(b)
        length(output) == n || throw(DimensionMismatch())
        for i = 1:n-1
            output[i+1] = b[i]
        end
        output[1] = b[n]
        output
    end

    # A' = cycle-forward operator
    function shiftfwd!(output, b)
        n = length(b)
        length(output) == n || throw(DimensionMismatch())
        for i = 2:n
            output[i-1] = b[i]
        end
        output[n] = b[1]
        output
    end

    b = rand(Float32, 5)
    Atimesb = [b[end]; b[1:end-1]]
    Aptimesb = [b[2:end]; b[1]]

    A = LinearMap(shiftback!, 5, 5, Int; ismutating=true)
    @fact eltype(A) --> Int
    @fact size(A) --> (5,5)
    @fact size(A,1) --> 5
    @fact size(A,2) --> 5
    @fact ndims(A) --> 2
    @fact length(A) --> 25
    @fact A*b --> Atimesb

    output = similar(b)
    @fact A_mul_B!(output, A, b) --> Atimesb
    @fact_throws A'*b

    A = LinearMap(shiftback!, shiftfwd!, 5, Int; ismutating=true)
    @fact eltype(A) --> Int
    @fact size(A) --> (5,5)
    @fact size(A,1) --> 5
    @fact size(A,2) --> 5
    @fact length(A) --> 25
    @fact A*b --> Atimesb

    output = similar(b)
    @fact A_mul_B!(output, A, b) --> Atimesb
    @fact A'*b --> Aptimesb
    @fact Ac_mul_B!(output, A, b) --> Aptimesb
end

end
