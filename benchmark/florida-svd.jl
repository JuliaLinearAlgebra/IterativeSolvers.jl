using IterativeSolvers
using MAT

BASEDIR = "florida"

for group in readdir(BASEDIR)
    isdir(joinpath(BASEDIR, group)) || continue
    for matrix in readdir(joinpath(BASEDIR, group))
        filename = joinpath(BASEDIR, group, matrix)
        print(filename, '\t')

        endswith(matrix, ".mat") || continue

        mf = matread(filename)
        if !haskey(mf, "Problem")
            warn("Skipping unknown file $filename: No 'Problem' struct")
            continue
        end

        pr = mf["Problem"]
        if !haskey(pr, "A")
            warn("Skipping unknown file $filename: 'Problem' struct has no matrix 'A'")
            continue
        end
        A = pr["A"]

        #eltype(A) <: Real || warn("Skipping matrix with unsupported element type $(eltype(A))")

        info("Size of matrix is $(size(A))")

        #info("Running naive SVD")
        #@time svdvals_gkl(A, 10)

        info("Thick restart")
        m, n = size(A)
        q = randn(n)
        eltype(A) <: Complex && (q += im*randn(n))
        scale!(q, inv(norm(q)))
        try
            @time svdvals_tr(A, q, 10)
        catch exc
            println("Exception: $exc")
        end
        #info("Harmonic restart")
        #@time svdvals_tr(A, q, 10, method=:harmonic)
        break
    end
end
