#Tests using the NIST Matrix Market

using IterativeSolvers, MatrixMarket
@info("Testing matrices from the NIST Matrix Market")
for (setname, matname) in [("cylshell", "s3dkq4m2")]
    dl_filename = string(setname, "_", filename)
    isfile(dl_filename) || download("ftp://math.nist.gov/pub/MatrixMarket2/misc/$setname/$matname.mtx.gz", dl_filename)

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
    #push!(solvers, gmres)
    for solver in solvers
        println("$dl_filename: Iterative solve using $solver")
        tic()
        xs, ch = solver(A, b)
        toc()
	r=norm(x-xs)
        println("Residual: $r matvecs: $(ch.mvps) iters: $(length(ch.residuals))")
    end
end

