#Tests using the University of Florida Sparse Matrix Collection
#http://www.cise.ufl.edu/research/sparse/matrices/

using IterativeSolvers, MAT

info("Testing matrices from the University of Florida Sparse Matrix Collection")
for (dirname, filename) in [("ACUSIM", "Pres_Poisson.mat")]
    dl_filename = string(dirname, "_", filename)
    isfile(dl_filename) || download("http://www.cise.ufl.edu/research/sparse/mat/$dirname/$filename", dl_filename)

    vars = matread(dl_filename)
    A, b = vars["Problem"]["A"], vars["Problem"]["b"][:]

    #Use direct solve as benchmark
    println("$dl_filename: Direct solve")
    tic()
    x = A \ b
    toc()

    @show Aissym = issym(A)
    @show Aisposdef = isposdef(A)
    solvers = {}
    Aissym && Aisposdef && push!(solvers, cg)
    push!(solvers, gmres)
    for solver in solvers
        println("$dl_filename: Iterative solve using $solver")
        tic()
        xs, ch = solver(A, b)
        toc()
	r=norm(x-xs)
        println("Residual: $r matvecs: $(ch.mvps) iters: $(length(ch.residuals))")
    end
end

