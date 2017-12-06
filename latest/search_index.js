var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#IterativeSolvers.jl-1",
    "page": "Home",
    "title": "IterativeSolvers.jl",
    "category": "section",
    "text": "IterativeSolvers.jl is a Julia package that provides efficient iterative algorithms for solving large linear systems, eigenproblems, and singular value problems. Most of the methods can be used matrix-free.For bug reports, feature requests and questions please submit an issue. If you're interested in contributing, please see the Contributing guide.For more information on future methods have a look at the package roadmap on deterministic methods, for randomized algorithms check here."
},

{
    "location": "index.html#What-method-should-I-use-for-linear-systems?-1",
    "page": "Home",
    "title": "What method should I use for linear systems?",
    "category": "section",
    "text": "When solving linear systems Ax = b for a square matrix A there are quite some options. The typical choices are listed below:Method When to use it\nConjugate Gradients Best choice for symmetric, positive-definite matrices\nMINRES For symmetric, indefinite matrices\nGMRES For nonsymmetric matrices when a good preconditioner is available\nIDR(s) For nonsymmetric, strongly indefinite problems without a good preconditioner\nBiCGStab(l) Otherwise for nonsymmetric problemsWe also offer Chebyshev iteration as an alternative to Conjugate Gradients when bounds on the spectrum are known.Stationary methods like Jacobi, Gauss-Seidel, SOR and SSOR can be used as smoothers to reduce high-frequency components in the error in just a few iterations.When solving least-squares problems we currently offer just LSMR and LSQR."
},

{
    "location": "index.html#Eigenproblems-and-SVD-1",
    "page": "Home",
    "title": "Eigenproblems and SVD",
    "category": "section",
    "text": "For the Singular Value Decomposition we offer SVDL, which is the Golub-Kahan-Lanczos procedure.For eigenvalue problems we have at this point just the Power Method and some convenience wrappers to do shift-and-invert."
},

{
    "location": "index.html#Randomized-algorithms-1",
    "page": "Home",
    "title": "Randomized algorithms",
    "category": "section",
    "text": "Randomized algorithms have gotten some traction lately. Some of the methods mentioned in [Halko2011] have been implemented as well, although their quality is generally poor. Also note that many classical methods such as the subspace iteration, BiCG and recent methods like IDR(s) are also \"randomized\" in some sense.[Halko2011]: Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. \"Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions.\" SIAM review 53.2 (2011): 217-288."
},

{
    "location": "getting_started.html#",
    "page": "Getting started",
    "title": "Getting started",
    "category": "page",
    "text": ""
},

{
    "location": "getting_started.html#Getting-started-1",
    "page": "Getting started",
    "title": "Getting started",
    "category": "section",
    "text": ""
},

{
    "location": "getting_started.html#Installation-1",
    "page": "Getting started",
    "title": "Installation",
    "category": "section",
    "text": "The package can be installed via Julia's package manager.julia> Pkg.add(\"IterativeSolvers\")"
},

{
    "location": "getting_started.html#Interface-1",
    "page": "Getting started",
    "title": "Interface",
    "category": "section",
    "text": "Virtually all solvers have the common function declarations:solver(A, args...; kwargs...)\nsolver!(x, A, args...; kwargs...)where A is a linear operator and x an initial guess. The second declaration also updates x in-place."
},

{
    "location": "getting_started.html#matrixfree-1",
    "page": "Getting started",
    "title": "Explicit matrices and the matrix-free approach",
    "category": "section",
    "text": "Rather than constructing an explicit matrix A of the type Matrix or SparseMatrixCSC, it is also possible to pass a general linear operator that performs matrix operations implicitly. This is called the matrix-free approach.For matrix-free types of A the following interface is expected to be defined:A*v computes the matrix-vector product on a v::AbstractVector;\nA_mul_B!(y, A, v) computes the matrix-vector product on a v::AbstractVector in-place;\neltype(A) returns the element type implicit in the equivalent matrix representation of A;\nsize(A, d) returns the nominal dimensions along the dth axis in the equivalent matrix representation of A.tip: Matrix-free with LinearMaps.jl\nWe strongly recommend LinearMaps.jl for matrix-free linear operators, as it implements the above methods already for you; you just have to write the action of the linear map."
},

{
    "location": "getting_started.html#Additional-arguments-1",
    "page": "Getting started",
    "title": "Additional arguments",
    "category": "section",
    "text": "Keyword names will vary depending on the method, however some of them will always have the same spelling:tol: (relative) stopping tolerance of the method;\nverbose: print information during the iterations;\nmaxiter: maximum number of allowed iterations;\nPl and Pr: left and right preconditioner. See Preconditioning;\nlog::Bool = false: output an extra element of type ConvergenceHistory containing the convergence history."
},

{
    "location": "getting_started.html#log-keyword-1",
    "page": "Getting started",
    "title": "log keyword",
    "category": "section",
    "text": "Most solvers contain the log keyword. This is to be used when obtaining more information is required, to use it place the set log to true.x, ch = cg(Master, rand(10, 10), rand(10) log=true)\nsvd, L, ch = svdl(Master, rand(100, 100), log=true)The function will now return one more parameter of type ConvergenceHistory."
},

{
    "location": "getting_started.html#IterativeSolvers.ConvergenceHistory",
    "page": "Getting started",
    "title": "IterativeSolvers.ConvergenceHistory",
    "category": "Type",
    "text": "Store general and in-depth information about an iterative method.\n\nFields\n\nmvps::Int: number of matrix vector products.\n\nmtvps::Int: number of transposed matrix-vector products\n\niters::Int: iterations taken by the method.\n\nrestart::T: restart relevant information.\n\nT == Int: iterations per restart.\nT == Void: methods without restarts.\n\nisconverged::Bool: convergence of the method.\n\ndata::Dict{Symbol,Any}: Stores all the information stored during the method execution. It stores tolerances, residuals and other information, e.g. Ritz values in svdl.\n\nConstructors\n\nConvergenceHistory()\nConvergenceHistory(restart)\n\nCreate ConvergenceHistory with empty fields.\n\nArguments\n\nrestart: number of iterations per restart.\n\nPlots\n\nSupports plots using the Plots.jl package via a type recipe. Vectors are ploted as series and matrices as scatterplots.\n\nImplements\n\nBase: getindex, setindex!, push!\n\n\n\n"
},

{
    "location": "getting_started.html#ConvergenceHistory-1",
    "page": "Getting started",
    "title": "ConvergenceHistory",
    "category": "section",
    "text": "A ConvergenceHistory instance stores information of a solver.Number of iterations.ch.itersConvergence status.ch.isconvergedStopping tolerances. (A Symbol key is needed to access)ch[:tol]Maximum number of iterations per restart. (Only on restarted methods)nrests(ch)Number of matrix-vectors and matrix-transposed-vector products.nprods(ch)Data stored on each iteration, accessed information can be either a vector or matrix. This data can be a lot of things, most commonly residual. (A Symbol key is needed to access)ch[:resnorm] #Vector or Matrix\nch[:resnorm, x] #Vector or Matrix element\nch[:resnorm, x, y] #Matrix elementConvergenceHistory"
},

{
    "location": "getting_started.html#Plotting-1",
    "page": "Getting started",
    "title": "Plotting",
    "category": "section",
    "text": "ConvergeHistory provides a recipe to use with the package Plots.jl, this makes it really easy to plot on different plot backends. There are two recipes provided:One for the whole ConvergenceHistory.plot(ch)The other one to plot data binded to a key._, ch = gmres(rand(10,10), rand(10), maxiter = 100, log=true)\nplot(ch, :resnorm, sep = :blue)Plot additional keywordssep::Symbol = :white: color of the line separator in restarted methods."
},

{
    "location": "preconditioning.html#",
    "page": "Preconditioning",
    "title": "Preconditioning",
    "category": "page",
    "text": ""
},

{
    "location": "preconditioning.html#Preconditioning-1",
    "page": "Preconditioning",
    "title": "Preconditioning",
    "category": "section",
    "text": "Many iterative solvers have the option to provide left and right preconditioners (Pl and Pr resp.) in order to speed up convergence or prevent stagnation. They transform a problem Ax = b into a better conditioned system (P_l^-1AP_r^-1)y = P_l^-1b, where x = P_r^-1y.These preconditioners should support the operations A_ldiv_B!(y, P, x) computes P \\ x in-place of y;\nA_ldiv_B!(P, x) computes P \\ x in-place of x;\nand P \\ x.If no preconditioners are passed to the solver, the method will default to Pl = Pr = IterativeSolvers.Identity()"
},

{
    "location": "preconditioning.html#Available-preconditioners-1",
    "page": "Preconditioning",
    "title": "Available preconditioners",
    "category": "section",
    "text": "IterativeSolvers.jl itself does not provide any other preconditioners besides Identity(), but recommends the following external packages:ILU.jl for incomplete LU decompositions (using drop tolerance);\nIncompleteSelectedInversion.jl for incomplete LDLt decompositions."
},

{
    "location": "linear_systems/cg.html#",
    "page": "Conjugate Gradients",
    "title": "Conjugate Gradients",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/cg.html#CG-1",
    "page": "Conjugate Gradients",
    "title": "Conjugate Gradients (CG)",
    "category": "section",
    "text": "Conjugate Gradients solves Ax = b approximately for x where A is a symmetric, positive-definite linear operator and b the right-hand side vector. The method uses short recurrences and therefore has fixed memory costs and fixed computational costs per iteration."
},

{
    "location": "linear_systems/cg.html#IterativeSolvers.cg",
    "page": "Conjugate Gradients",
    "title": "IterativeSolvers.cg",
    "category": "Function",
    "text": "cg(A, b; kwargs...) -> x, [history]\n\nSame as cg!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/cg.html#IterativeSolvers.cg!",
    "page": "Conjugate Gradients",
    "title": "IterativeSolvers.cg!",
    "category": "Function",
    "text": "cg!(x, A, b; kwargs...) -> x, [history]\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool: If true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\nPl = Identity(): left preconditioner of the method. Should be symmetric,  positive-definite like A.\ntol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol;\nmaxiter::Int = size(A,2): maximum number of iterations;\nverbose::Bool = false: print method information;\nlog::Bool = false: keep track of the residual norm in each iteration;\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "linear_systems/cg.html#Usage-1",
    "page": "Conjugate Gradients",
    "title": "Usage",
    "category": "section",
    "text": "cg\ncg!"
},

{
    "location": "linear_systems/cg.html#Implementation-details-1",
    "page": "Conjugate Gradients",
    "title": "Implementation details",
    "category": "section",
    "text": "The current implementation follows a rather standard approach. Note that preconditioned CG (or PCG) is slightly different from ordinary CG, because the former must compute the residual explicitly, while it is available as byproduct in the latter. Our implementation of CG ensures the minimal number of vector operations.tip: Tip\nCG can be used as an iterator."
},

{
    "location": "linear_systems/chebyshev.html#",
    "page": "Chebyshev iteration",
    "title": "Chebyshev iteration",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/chebyshev.html#Chebyshev-1",
    "page": "Chebyshev iteration",
    "title": "Chebyshev iteration",
    "category": "section",
    "text": "Chebyshev iteration solves the problem Ax=b approximately for x where A is a symmetric, definite linear operator and b the right-hand side vector. The methods assumes the interval lambda_min lambda_max containing all eigenvalues of A is known, so that x can be iteratively constructed via a Chebyshev polynomial with zeros in this interval. This polynomial ultimately acts as a filter that removes components in the direction of the eigenvectors from the initial residual.The main advantage with respect to Conjugate Gradients is that BLAS1 operations such as inner products are avoided."
},

{
    "location": "linear_systems/chebyshev.html#IterativeSolvers.chebyshev",
    "page": "Chebyshev iteration",
    "title": "IterativeSolvers.chebyshev",
    "category": "Function",
    "text": "chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]\n\nSame as chebyshev!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/chebyshev.html#IterativeSolvers.chebyshev!",
    "page": "Chebyshev iteration",
    "title": "IterativeSolvers.chebyshev!",
    "category": "Function",
    "text": "chebyshev!(x, A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]\n\nSolve Ax = b for symmetric, definite matrices A using Chebyshev iteration.\n\nArguments\n\nx: initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side;\nλmin::Real: lower bound for the real eigenvalues\nλmax::Real: upper bound for the real eigenvalues\n\nKeywords\n\ninitially_zero::Bool = false: if true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\ntol: tolerance for stopping condition |r_k| / |r_0| ≤ tol.\nmaxiter::Int = size(A, 2): maximum number of inner iterations of GMRES;\nPl = Identity(): left preconditioner;\nlog::Bool = false: keep track of the residual norm in each iteration;\nverbose::Bool = false: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "linear_systems/chebyshev.html#Usage-1",
    "page": "Chebyshev iteration",
    "title": "Usage",
    "category": "section",
    "text": "chebyshev\nchebyshev!"
},

{
    "location": "linear_systems/chebyshev.html#Implementation-details-1",
    "page": "Chebyshev iteration",
    "title": "Implementation details",
    "category": "section",
    "text": "warning: BLAS1 operations\nAlthough the method is often used to avoid computation of inner products, the stopping criterion is still based on the residual norm. Hence the current implementation is not free of BLAS1 operations.tip: Tip\nChebyshev iteration can be used as an iterator."
},

{
    "location": "linear_systems/minres.html#",
    "page": "MINRES",
    "title": "MINRES",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/minres.html#MINRES-1",
    "page": "MINRES",
    "title": "MINRES",
    "category": "section",
    "text": "MINRES is a short-recurrence version of GMRES for solving Ax = b approximately for x where A is a symmetric, Hermitian, skew-symmetric or skew-Hermitian linear operator and b the right-hand side vector."
},

{
    "location": "linear_systems/minres.html#IterativeSolvers.minres",
    "page": "MINRES",
    "title": "IterativeSolvers.minres",
    "category": "Function",
    "text": "minres(A, b; kwargs...) -> x, [history]\n\nSame as minres!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/minres.html#IterativeSolvers.minres!",
    "page": "MINRES",
    "title": "IterativeSolvers.minres!",
    "category": "Function",
    "text": "minres!(x, A, b; kwargs...) -> x, [history]\n\nSolve Ax = b for (skew-)Hermitian matrices A using MINRES.\n\nArguments\n\nx: initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool = false: if true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\nskew_hermitian::Bool = false: if true assumes that A is skew-symmetric or skew-Hermitian;\ntol: tolerance for stopping condition |r_k| / |r_0| ≤ tol. Note that the residual is computed only approximately;\nmaxiter::Int = size(A, 2): maximum number of iterations;\nPl: left preconditioner;\nPr: right preconditioner;\nlog::Bool = false: keep track of the residual norm in each iteration;\nverbose::Bool = false: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "linear_systems/minres.html#Usage-1",
    "page": "MINRES",
    "title": "Usage",
    "category": "section",
    "text": "minres\nminres!"
},

{
    "location": "linear_systems/minres.html#Implementation-details-1",
    "page": "MINRES",
    "title": "Implementation details",
    "category": "section",
    "text": "MINRES exploits the tridiagonal structure of the Hessenberg matrix. Although MINRES is mathematically equivalent to GMRES, it might not be equivalent in finite precision. MINRES updates the solution asx = x_0 + (V R^-1) (Q^*r_0e_1)where V is the orthonormal basis for the Krylov subspace and QR is the QR-decomposition of the Hessenberg matrix. Note that the brackets are placed slightly differently from how GMRES would update the residual.MINRES computes V and W = VR^-1 via a three-term recurrence, using only the last column of R Therefore we pre-allocate only six vectors, save only the last two entries of Q^*r_0e_1 and part of the last column of the Hessenberg matrix.note: Real and complex arithmetic\nIf A is Hermitian, then the Hessenberg matrix will be real. This is exploited in the current implementation.If A is skew-Hermitian, the diagonal of the Hessenberg matrix will be imaginary, and hence we use complex arithmetic in that case.tip: Tip\nMINRES can be used as an iterator."
},

{
    "location": "linear_systems/bicgstabl.html#",
    "page": "BiCGStab(l)",
    "title": "BiCGStab(l)",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/bicgstabl.html#BiCGStabl-1",
    "page": "BiCGStab(l)",
    "title": "BiCGStab(l)",
    "category": "section",
    "text": "BiCGStab(l) solves the problem Ax = b approximately for x where A is a general, linear operator and b the right-hand side vector. The methods combines BiCG with l GMRES iterations, resulting in a short-reccurence iteration. As a result the memory is fixed as well as the computational costs per iteration."
},

{
    "location": "linear_systems/bicgstabl.html#IterativeSolvers.bicgstabl",
    "page": "BiCGStab(l)",
    "title": "IterativeSolvers.bicgstabl",
    "category": "Function",
    "text": "bicgstabl(A, b, l; kwargs...) -> x, [history]\n\nSame as bicgstabl!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/bicgstabl.html#IterativeSolvers.bicgstabl!",
    "page": "BiCGStab(l)",
    "title": "IterativeSolvers.bicgstabl!",
    "category": "Function",
    "text": "bicgstabl!(x, A, b, l; kwargs...) -> x, [history]\n\nArguments\n\nA: linear operator;\nb: right hand side (vector);\nl::Int = 2: Number of GMRES steps.\n\nKeywords\n\nmax_mv_products::Int = size(A, 2): maximum number of matrix vector products.\n\nFor BiCGStab(l) this is a less dubious term than \"number of iterations\";\n\nPl = Identity(): left preconditioner of the method;\ntol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol.   Note that (1) the true residual norm is never computed during the iterations,   only an approximation; and (2) if a preconditioner is given, the stopping condition is based on the   preconditioned residual.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "linear_systems/bicgstabl.html#Usage-1",
    "page": "BiCGStab(l)",
    "title": "Usage",
    "category": "section",
    "text": "bicgstabl\nbicgstabl!"
},

{
    "location": "linear_systems/bicgstabl.html#Implementation-details-1",
    "page": "BiCGStab(l)",
    "title": "Implementation details",
    "category": "section",
    "text": "The method is based on the original article [Sleijpen1993], but does not implement later improvements. The normal equations arising from the GMRES steps are solved without orthogonalization. Hence the method should only be reliable for relatively small values of l.The r and u factors are pre-allocated as matrices of size n times (l + 1), so that BLAS2 methods can be used. Also the random shadow residual is pre-allocated as a vector. Hence the storage costs are approximately 2l + 3 vectors.tip: Tip\nBiCGStabl(l) can be used as an iterator.[Sleijpen1993]: Sleijpen, Gerard LG, and Diederik R. Fokkema. \"BiCGstab(l) for  linear equations involving unsymmetric matrices with complex spectrum.\"  Electronic Transactions on Numerical Analysis 1.11 (1993): 2000."
},

{
    "location": "linear_systems/idrs.html#",
    "page": "IDR(s)",
    "title": "IDR(s)",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/idrs.html#IDRs-1",
    "page": "IDR(s)",
    "title": "IDR(s)",
    "category": "section",
    "text": "The Induced Dimension Reduction method is a family of simple and fast Krylov subspace algorithms for solving large nonsymmetric linear systems. The idea behind the IDR(s) variant is to generate residuals that are in the nested subspaces of shrinking dimensions."
},

{
    "location": "linear_systems/idrs.html#IterativeSolvers.idrs",
    "page": "IDR(s)",
    "title": "IterativeSolvers.idrs",
    "category": "Function",
    "text": "idrs(A, b; s = 8) -> x, [history]\n\nSame as idrs!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/idrs.html#IterativeSolvers.idrs!",
    "page": "IDR(s)",
    "title": "IterativeSolvers.idrs!",
    "category": "Function",
    "text": "idrs!(x, A, b; s = 8) -> x, [history]\n\nSolve the problem Ax = b approximately with IDR(s), where s is the dimension of the shadow space.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ns::Integer = 8: dimension of the shadow space;\ntol: relative tolerance;\nmaxiter::Int = size(A, 2): maximum number of iterations;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "linear_systems/idrs.html#Usage-1",
    "page": "IDR(s)",
    "title": "Usage",
    "category": "section",
    "text": "idrs\nidrs!"
},

{
    "location": "linear_systems/idrs.html#Implementation-details-1",
    "page": "IDR(s)",
    "title": "Implementation details",
    "category": "section",
    "text": "The current implementation is based on the MATLAB version by Van Gijzen and Sonneveld. For background see [Sonneveld2008], [VanGijzen2011] and the IDR(s) webpage.[Sonneveld2008]: IDR(s): a family of simple and fast algorithms for solving large nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035–1062, 2008[VanGijzen2011]: Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011"
},

{
    "location": "linear_systems/gmres.html#",
    "page": "Restarted GMRES",
    "title": "Restarted GMRES",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/gmres.html#GMRES-1",
    "page": "Restarted GMRES",
    "title": "Restarted GMRES",
    "category": "section",
    "text": "GMRES solves the problem Ax = b approximately for x where A is a general, linear operator and b the right-hand side vector. The method is optimal in the sense that it selects the solution with minimal residual from a Krylov subspace, but the price of optimality is increasing storage and computation effort per iteration. Restarts are necessary to fix these costs."
},

{
    "location": "linear_systems/gmres.html#IterativeSolvers.gmres",
    "page": "Restarted GMRES",
    "title": "IterativeSolvers.gmres",
    "category": "Function",
    "text": "gmres(A, b; kwargs...) -> x, [history]\n\nSame as gmres!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/gmres.html#IterativeSolvers.gmres!",
    "page": "Restarted GMRES",
    "title": "IterativeSolvers.gmres!",
    "category": "Function",
    "text": "gmres!(x, A, b; kwargs...) -> x, [history]\n\nSolves the problem Ax = b with restarted GMRES.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool: If true assumes that iszero(x) so that one  matrix-vector product can be saved when computing the initial  residual vector;\ntol: relative tolerance;\nrestart::Int = min(20, size(A, 2)): restarts GMRES after specified number of iterations;\nmaxiter::Int = size(A, 2): maximum number of inner iterations of GMRES;\nPl: left preconditioner;\nPr: right preconditioner;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n"
},

{
    "location": "linear_systems/gmres.html#Usage-1",
    "page": "Restarted GMRES",
    "title": "Usage",
    "category": "section",
    "text": "gmres\ngmres!"
},

{
    "location": "linear_systems/gmres.html#Implementation-details-1",
    "page": "Restarted GMRES",
    "title": "Implementation details",
    "category": "section",
    "text": "The implementation pre-allocates a matrix V of size n by restart whose columns form an orthonormal basis for the Krylov subspace. This allows BLAS2 operations when updating the solution vector x. The Hessenberg matrix is also pre-allocated.Modified Gram-Schmidt is used to orthogonalize the columns of V.The computation of the residual norm is implemented in a non-standard way, namely keeping track of a vector gamma in the null-space of H_k^*, which is the adjoint of the (k + 1) times k Hessenberg matrix H_k at the kth iteration. Only when x needs to be updated is the Hessenberg matrix mutated with Givens rotations.tip: Tip\nGMRES can be used as an iterator. This makes it possible to access the Hessenberg matrix and Krylov basis vectors during the iterations."
},

{
    "location": "linear_systems/lsmr.html#",
    "page": "LSMR",
    "title": "LSMR",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/lsmr.html#LSMR-1",
    "page": "LSMR",
    "title": "LSMR",
    "category": "section",
    "text": "Least-squares minimal residual"
},

{
    "location": "linear_systems/lsmr.html#IterativeSolvers.lsmr",
    "page": "LSMR",
    "title": "IterativeSolvers.lsmr",
    "category": "Function",
    "text": "lsmr(A, b; kwrags...) -> x, [history]\n\nSame as lsmr!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/lsmr.html#IterativeSolvers.lsmr!",
    "page": "LSMR",
    "title": "IterativeSolvers.lsmr!",
    "category": "Function",
    "text": "lsmr!(x, A, b; kwargs...) -> x, [history]\n\nMinimizes Ax - b^2 + x^2 in the Euclidean norm. If multiple solutions exists the minimum norm solution is returned.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is  algebraically equivalent to applying MINRES to the normal equations  (A^*A + ^2I)x = A^*b, but has better numerical properties,  especially if A is ill-conditioned.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\nconlim::Number = 1e8: stopping tolerance. lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting\natol = btol = conlim = zero, but the number of iterations may then be excessive.\nmaxiter::Int = maximum(size(A)): maximum number of iterations.\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n:btol => ::Real: btol stopping tolerance.\n:ctol => ::Real: ctol stopping tolerance.\n:anorm => ::Real: anorm.\n:rnorm => ::Real: rnorm.\n:cnorm => ::Real: cnorm.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "linear_systems/lsmr.html#Usage-1",
    "page": "LSMR",
    "title": "Usage",
    "category": "section",
    "text": "lsmr\nlsmr!"
},

{
    "location": "linear_systems/lsmr.html#Implementation-details-1",
    "page": "LSMR",
    "title": "Implementation details",
    "category": "section",
    "text": "Adapted from: http://web.stanford.edu/group/SOL/software/lsmr/"
},

{
    "location": "linear_systems/lsqr.html#",
    "page": "LSQR",
    "title": "LSQR",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/lsqr.html#LSQR-1",
    "page": "LSQR",
    "title": "LSQR",
    "category": "section",
    "text": ""
},

{
    "location": "linear_systems/lsqr.html#IterativeSolvers.lsqr",
    "page": "LSQR",
    "title": "IterativeSolvers.lsqr",
    "category": "Function",
    "text": "lsqr(A, b; kwrags...) -> x, [history]\n\nSame as lsqr!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/lsqr.html#IterativeSolvers.lsqr!",
    "page": "LSQR",
    "title": "IterativeSolvers.lsqr!",
    "category": "Function",
    "text": "lsqr!(x, A, b; kwargs...) -> x, [history]\n\nMinimizes Ax - b^2 + damp*x^2 in the Euclidean norm. If multiple solutions exists returns the minimal norm solution.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying CG to the normal equations  (A^*A + ^2I)x = A^*b but has better numerical properties,  especially if A is ill-conditioned.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ndamp::Number = 0: damping parameter.\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\nmaxiter::Int = maximum(size(A)): maximum number of iterations.\nverbose::Bool = false: print method information.\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nReturn values\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n:btol => ::Real: btol stopping tolerance.\n:ctol => ::Real: ctol stopping tolerance.\n:anorm => ::Real: anorm.\n:rnorm => ::Real: rnorm.\n:cnorm => ::Real: cnorm.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "linear_systems/lsqr.html#Usage-1",
    "page": "LSQR",
    "title": "Usage",
    "category": "section",
    "text": "lsqr\nlsqr!"
},

{
    "location": "linear_systems/lsqr.html#Implementation-details-1",
    "page": "LSQR",
    "title": "Implementation details",
    "category": "section",
    "text": "Adapted from: http://web.stanford.edu/group/SOL/software/lsqr/."
},

{
    "location": "linear_systems/stationary.html#",
    "page": "Stationary methods",
    "title": "Stationary methods",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/stationary.html#Stationary-1",
    "page": "Stationary methods",
    "title": "Stationary methods",
    "category": "section",
    "text": "Stationary methods are typically used as smoothers in multigrid methods, where only very few iterations are applied to get rid of high-frequency components in the error. The implementations of stationary methods have this goal in mind, which means there is no other stopping criterion besides the maximum number of iterations.note: CSC versus CSR\nJulia stores matrices column-major. In order to avoid cache misses, the implementations of our stationary methods traverse the matrices column-major. This deviates from classical textbook implementations. Also the SOR and SSOR methods cannot be computed efficiently in-place, but require a temporary vector.When it comes to SparseMatrixCSC, we precompute in all stationary methods an integer array of the indices of the diagonal to avoid expensive searches in each iteration."
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.jacobi",
    "page": "Stationary methods",
    "title": "IterativeSolvers.jacobi",
    "category": "Function",
    "text": "jacobi(A, b) -> x\n\nSame as jacobi!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.jacobi!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.jacobi!",
    "category": "Function",
    "text": "jacobi!(x, A::AbstractMatrix, b; maxiter=10) -> x\n\nPerforms exactly maxiter Jacobi iterations.\n\nAllocates a single temporary vector and traverses A columnwise.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\njacobi!(x, A::SparseMatrixCSC, b; maxiter=10) -> x\n\nPerforms exactly maxiter Jacobi iterations.\n\nAllocates a temporary vector and precomputes the diagonal indices.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#Jacobi-1",
    "page": "Stationary methods",
    "title": "Jacobi",
    "category": "section",
    "text": "jacobi\njacobi!"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.gauss_seidel",
    "page": "Stationary methods",
    "title": "IterativeSolvers.gauss_seidel",
    "category": "Function",
    "text": "gauss_seidel(A, b) -> x\n\nSame as gauss_seidel!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.gauss_seidel!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.gauss_seidel!",
    "category": "Function",
    "text": "gauss_seidel!(x, A::AbstractMatrix, b; maxiter=10) -> x\n\nPerforms exactly maxiter Gauss-Seidel iterations.\n\nWorks fully in-place and traverses A columnwise.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\ngauss_seidel!(x, A::SparseMatrixCSC, b; maxiter=10) -> x\n\nPerforms exactly maxiter Gauss-Seidel iterations.\n\nWorks fully in-place, but precomputes the diagonal indices.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#Gauss-Seidel-1",
    "page": "Stationary methods",
    "title": "Gauss-Seidel",
    "category": "section",
    "text": "gauss_seidel\ngauss_seidel!"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.sor",
    "page": "Stationary methods",
    "title": "IterativeSolvers.sor",
    "category": "Function",
    "text": "sor(A, b, ω::Real) -> x\n\nSame as sor!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.sor!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.sor!",
    "category": "Function",
    "text": "sor!(x, A::AbstractMatrix, b, ω::Real; maxiter=10) -> x\n\nPerforms exactly maxiter SOR iterations with relaxation parameter ω.\n\nAllocates a single temporary vector and traverses A columnwise.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\nsor!(x, A::SparseMatrixCSC, b, ω::Real; maxiter=10)\n\nPerforms exactly maxiter SOR iterations with relaxation parameter ω.\n\nAllocates a temporary vector and precomputes the diagonal indices.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#SOR-1",
    "page": "Stationary methods",
    "title": "Successive over-relaxation (SOR)",
    "category": "section",
    "text": "sor\nsor!"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.ssor",
    "page": "Stationary methods",
    "title": "IterativeSolvers.ssor",
    "category": "Function",
    "text": "ssor(A, b, ω::Real) -> x\n\nSame as ssor!, but allocates a solution vector x initialized with zeros.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#IterativeSolvers.ssor!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.ssor!",
    "category": "Function",
    "text": "ssor!(x, A::AbstractMatrix, b, ω::Real; maxiter=10) -> x\n\nPerforms exactly maxiter SSOR iterations with relaxation parameter ω. Each iteration  is basically a forward and backward sweep of SOR.\n\nAllocates a single temporary vector and traverses A columnwise.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\nssor!(x, A::SparseMatrixCSC, b, ω::Real; maxiter=10)\n\nPerforms exactly maxiter SSOR iterations with relaxation parameter ω. Each iteration  is basically a forward and backward sweep of SOR.\n\nAllocates a temporary vector and precomputes the diagonal indices.\n\nThrows Base.LinAlg.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n"
},

{
    "location": "linear_systems/stationary.html#SSOR-1",
    "page": "Stationary methods",
    "title": "Symmetric successive over-relaxation (SSOR)",
    "category": "section",
    "text": "ssor\nssor!tip: Tip\nAll stationary methods can be used a iterators."
},

{
    "location": "eigenproblems/power_method.html#",
    "page": "Power method",
    "title": "Power method",
    "category": "page",
    "text": ""
},

{
    "location": "eigenproblems/power_method.html#PowerMethod-1",
    "page": "Power method",
    "title": "(Inverse) power method",
    "category": "section",
    "text": "Solves the eigenproblem Ax = x approximately where A is a general linear map. By default converges towards the dominant eigenpair ( x) such that  is largest. Shift-and-invert can be applied to target a specific eigenvalue near shift in the complex plane."
},

{
    "location": "eigenproblems/power_method.html#IterativeSolvers.powm",
    "page": "Power method",
    "title": "IterativeSolvers.powm",
    "category": "Function",
    "text": "powm(B; kwargs...) -> λ, x, [history]\n\nSee powm!. Calls powm!(B, x0; kwargs...) with  x0 initialized as a random, complex unit vector.\n\n\n\n"
},

{
    "location": "eigenproblems/power_method.html#IterativeSolvers.powm!",
    "page": "Power method",
    "title": "IterativeSolvers.powm!",
    "category": "Function",
    "text": "powm!(B, x; shift = zero(eltype(B)), inverse::Bool = false, kwargs...) -> λ, x, [history]\n\nBy default finds the approximate eigenpair (λ, x) of B where |λ| is largest.\n\nArguments\n\nB: linear map, see the note below.\nx: normalized initial guess. Don't forget to use complex arithmetic when necessary.\n\nKeywords\n\ntol::Real = eps(real(eltype(B))) * size(B, 2) ^ 3: stopping tolerance for the residual norm;\nmaxiter::Integer = size(B,2): maximum number of iterations;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nnote: Shift-and-invert\nWhen applying shift-and-invert to Ax = x with invert = true and shift = ..., note  that the role of B * b becomes computing inv(A - shift I) * b. So rather than  passing the linear map A itself, pass a linear map B that has the action of  shift-and-invert. The eigenvalue is transformed back to an eigenvalue of the actual  matrix A.\n\nReturn values\n\nif log is false\n\nλ::Number approximate eigenvalue computed as the Rayleigh quotient;\nx::Vector approximate eigenvector.\n\nif log is true\n\nλ::Number: approximate eigenvalue computed as the Rayleigh quotient;\nx::Vector: approximate eigenvector;\nhistory: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance;\n:resnom => ::Vector: residual norm at each iteration.\n\nExamples\n\nusing LinearMaps\nσ = 1.0 + 1.3im\nA = rand(Complex128, 50, 50)\nF = lufact(A - σ * I)\nFmap = LinearMap{Complex128}((y, x) -> A_ldiv_B!(y, F, x), 50, ismutating = true)\nλ, x = powm(Fmap, inverse = true, shift = σ, tol = 1e-4, maxiter = 200)\n\n\n\n"
},

{
    "location": "eigenproblems/power_method.html#IterativeSolvers.invpowm",
    "page": "Power method",
    "title": "IterativeSolvers.invpowm",
    "category": "Function",
    "text": "invpowm(B; shift = σ, kwargs...) -> λ, x, [history]\n\nFind the approximate eigenpair (λ, x) of A near shift, where B is a linear map that has the effect B * v = inv(A - σI) * v.\n\nThe method calls powm!(B, x0; inverse = true, shift = σ) with x0 a random, complex unit vector. See powm!\n\nExamples\n\nusing LinearMaps\nσ = 1.0 + 1.3im\nA = rand(Complex128, 50, 50)\nF = lufact(A - σ * I)\nFmap = LinearMap{Complex128}((y, x) -> A_ldiv_B!(y, F, x), 50, ismutating = true)\nλ, x = invpowm(Fmap, shift = σ, tol = 1e-4, maxiter = 200)\n\n\n\n"
},

{
    "location": "eigenproblems/power_method.html#IterativeSolvers.invpowm!",
    "page": "Power method",
    "title": "IterativeSolvers.invpowm!",
    "category": "Function",
    "text": "invpowm!(B, x0; shift = σ, kwargs...) -> λ, x, [history]\n\nFind the approximate eigenpair (λ, x) of A near shift, where B is a linear map that has the effect B * v = inv(A - σI) * v.\n\nThe method calls powm!(B, x0; inverse = true, shift = σ). See powm!.\n\n\n\n"
},

{
    "location": "eigenproblems/power_method.html#Usage-1",
    "page": "Power method",
    "title": "Usage",
    "category": "section",
    "text": "powm\npowm!\ninvpowm\ninvpowm!"
},

{
    "location": "eigenproblems/power_method.html#Implementation-details-1",
    "page": "Power method",
    "title": "Implementation details",
    "category": "section",
    "text": "Storage requirements are 3 vectors: the approximate eigenvector x, the residual vector r and a temporary. The residual norm lags behind one iteration, as it is computed when Ax is performed. Therefore the final resdiual norm is even smaller."
},

{
    "location": "svd/svdl.html#",
    "page": "SVDL",
    "title": "SVDL",
    "category": "page",
    "text": ""
},

{
    "location": "svd/svdl.html#SVDL-1",
    "page": "SVDL",
    "title": "Golub-Kahan-Lanczos (SVDL)",
    "category": "section",
    "text": "The SVDL method computes a partial, approximate SVD decomposition of a general linear operator A."
},

{
    "location": "svd/svdl.html#IterativeSolvers.svdl",
    "page": "SVDL",
    "title": "IterativeSolvers.svdl",
    "category": "Function",
    "text": "svdl(A) -> Σ, L, [history]\n\nCompute some singular values (and optionally vectors) using Golub-Kahan-Lanczos bidiagonalization [Golub1965] with thick restarting [Wu2000].\n\nIf log is set to true is given, method will output a tuple X, L, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return X, L.\n\nArguments\n\nA : The matrix or matrix-like object whose singular values are desired.\n\nKeywords\n\nnsv::Int = 6: number of singular values requested;\nv0 = random unit vector: starting guess vector in the domain of A. The length of q should be the number of columns in A;\nk::Int = 2nsv: maximum number of Lanczos vectors to compute before restarting;\nj::Int = nsv: number of vectors to keep at the end of the restart.  We don't recommend j < nsv;\nmaxiter::Int = minimum(size(A)): maximum number of iterations to run;\nverbose::Bool = false: print information at each iteration;\ntol::Real = √eps(): maximum absolute error in each desired singular value;\nreltol::Real=√eps(): maximum error in each desired singular value relative to the  estimated norm of the input matrix;\nmethod::Symbol=:ritz: restarting algorithm to use. Valid choices are:\n:ritz: Thick restart with Ritz values [Wu2000].\n:harmonic: Restart with harmonic Ritz values [Baglama2005].\nvecs::Symbol = :none: singular vectors to return.\n:both: Both left and right singular vectors are returned.\n:left: Only the left singular vectors are returned.\n:right: Only the right singular vectors are returned.\n:none: No singular vectors are returned.\ndolock::Bool=false: If true, locks converged Ritz values, removing them from the Krylov subspace being searched in the next macroiteration;\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nReturn values\n\nif log is false\n\nΣ: list of the desired singular values if vecs == :none (the default), otherwise  returns an SVD object with the desired singular vectors filled in;\nL: computed partial factorizations of A.\n\nif log is true\n\nΣ: list of the desired singular values if vecs == :none (the default),\n\notherwise returns an SVD object with the desired singular vectors filled in;\n\nL: computed partial factorizations of A;\nhistory: convergence history.\n\nConvergenceHistory keys\n\n:betas => betas: The history of the computed betas.\n:Bs => Bs: The history of the computed projected matrices.\n:ritz => ritzvalhist: Ritz values computed at each iteration.\n:conv => convhist: Convergence data.\n\n[Golub1965]: Golub, Gene, and William Kahan. \"Calculating the singular values and pseudo-inverse  of a matrix.\" Journal of the Society for Industrial and Applied Mathematics,  Series B: Numerical Analysis 2.2 (1965): 205-224.\n\n[Wu2000]: Wu, Kesheng, and Horst Simon. \"Thick-restart Lanczos method for large symmetric  eigenvalue problems.\" SIAM Journal on Matrix Analysis and Applications 22.2  (2000): 602-616.\n\n\n\n"
},

{
    "location": "svd/svdl.html#Usage-1",
    "page": "SVDL",
    "title": "Usage",
    "category": "section",
    "text": "svdl"
},

{
    "location": "svd/svdl.html#Implementation-details-1",
    "page": "SVDL",
    "title": "Implementation details",
    "category": "section",
    "text": "The implementation of thick restarting follows closely that of SLEPc as described in [Hernandez2008]. Thick restarting can be turned off by setting k = maxiter, but most of the time this is not desirable.The singular vectors are computed directly by forming the Ritz vectors from the product of the Lanczos vectors L.P/L.Q and the singular vectors of L.B. Additional accuracy in the singular triples can be obtained using inverse iteration.[Hernandez2008]: Vicente Hernández, José E. Román, and Andrés Tomás. \"A robust and efficient parallel SVD solver based on restarted Lanczos bidiagonalization.\" Electronic Transactions on Numerical Analysis 31 (2008): 68-85."
},

{
    "location": "randomized.html#",
    "page": "Randomized algorithms",
    "title": "Randomized algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "randomized.html#IterativeSolvers.reig",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.reig",
    "category": "Function",
    "text": "reig(A, l)\n\nCompute the spectral (Eigen) decomposition of A using a randomized algorithm.\n\nArguments\n\nA: input matrix.\nl::Int: number of eigenpairs to find.\n\nOutput\n\n::Base.LinAlg.Eigen: eigen decomposition.\n\nImplementation note\n\nThis is a wrapper around eigfact_onepass() which uses the randomized samples found using srft(l).\n\n\n\n"
},

{
    "location": "randomized.html#IterativeSolvers.rsvdfact",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.rsvdfact",
    "category": "Function",
    "text": "rsvdfact(A, n, p=0)\n\nCompute partial singular value decomposition of A using a randomized algorithm.\n\nArguments\n\nA: input matrix.\n\nn::Int: number of singular value/vector pairs to find.\n\np::Int=0: number of extra vectors to include in computation.\n\nOutput\n\n::SVD: singular value decomposition.\n\nwarning: Accuracy\nThis variant of the randomized singular value decomposition is the most commonly found implementation but is not recommended for accurate computations, as it often has trouble finding the n largest singular pairs, but rather finds n large singular pairs which may not necessarily be the largest.\n\nImplementation note\n\nThis function calls rrange, which uses naive randomized rangefinding to compute a basis for a subspace of dimension n (Algorithm 4.1 of [Halko2011]), followed by svdfact_restricted(), which computes the exact SVD factorization on the restriction of A to this randomly selected subspace (Algorithm 5.1 of [Halko2011]).\n\nAlternatively, you can mix and match your own randomized algorithm using any of the randomized range finding algorithms to find a suitable subspace and feeding the result to one of the routines that computes the SVD restricted to that subspace.\n\n\n\n"
},

{
    "location": "randomized.html#IterativeSolvers.rsvd_fnkz",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.rsvd_fnkz",
    "category": "Function",
    "text": "rsvd_fnkz(A, k)\n\nCompute the randomized SVD by iterative refinement from randomly selected columns/rows [Friedland2006].\n\nArguments\n\nA: matrix whose SVD is desired;\nk::Int: desired rank of approximation (k ≤ min(m, n)).\n\nKeywords\n\nl::Int = k: number of columns/rows to sample at each iteration (1 ≤ l ≤ k);\nN::Int = minimum(size(A)): maximum number of iterations;\nϵ::Real = prod(size(A))*eps(): relative threshold for convergence, as measured by growth of the spectral norm;\nmethod::Symbol = :eig: problem to solve.\n:eig: eigenproblem.\n:svd: singular problem.\nverbose::Bool = false: print convergence information at each iteration.\n\nReturn value\n\nSVD object of rank ≤ k.\n\n[Friedland2006]: Friedland, Shmuel, et al. \"Fast Monte-Carlo low rank approximations for  matrices.\" System of Systems Engineering, 2006 IEEE/SMC International  Conference on. IEEE, 2006.\n\n\n\n"
},

{
    "location": "randomized.html#Randomized-1",
    "page": "Randomized algorithms",
    "title": "Randomized algorithms",
    "category": "section",
    "text": "The methods below are based on [Halko2011].reig\nrsvdfact\nrsvd_fnkz"
},

{
    "location": "randomized.html#IterativeSolvers.rcond",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.rcond",
    "category": "Function",
    "text": "rcond(A, iters=1)\n\nEstimate matrix condition number randomly.\n\nArguments\n\nA: matrix whose condition number to estimate. Must be square and\n\nsupport premultiply (A*⋅) and solve (A\\⋅).\n\niters::Int = 1: number of power iterations to run.\n\nKeywords\n\np::Real = 0.05: probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains κ(A) with probability 1 - p.\n\nImplementation note\n\n[Dixon1983] originally describes this as a computation that can be done by computing the necessary number of power iterations given p and the desired accuracy parameter θ=y/x. However, these bounds were only derived under the assumptions of exact arithmetic. Empirically, iters≥4 has been seen to result in incorrect results in that the computed interval does not contain the true condition number. This implemention therefore makes iters an explicitly user-controllable parameter from which to infer the accuracy parameter and hence the interval containing κ(A). ```\n\n\n\n"
},

{
    "location": "randomized.html#Condition-number-estimate-1",
    "page": "Randomized algorithms",
    "title": "Condition number estimate",
    "category": "section",
    "text": "rcond"
},

{
    "location": "randomized.html#IterativeSolvers.reigmin",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.reigmin",
    "category": "Function",
    "text": "reigmin(A, iters=1)\n\nEstimate minimal eigenvalue randomly.\n\nArguments\n\nA: Matrix whose maximal eigenvalue to estimate.\n\nMust be square and support premultiply (A*⋅).\n\niters::Int=1: Number of power iterations to run. (Recommended: iters ≤ 3)\n\nKeywords\n\np::Real=0.05: Probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains the maximal eigenvalue of A with probability 1 - p.\n\nReferences\n\nCorollary of Theorem 1 in [Dixon1983]\n\n\n\n"
},

{
    "location": "randomized.html#IterativeSolvers.reigmax",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.reigmax",
    "category": "Function",
    "text": "reigmax(A, iters=1)\n\nEstimate maximal eigenvalue randomly.\n\nArguments\n\nA: Matrix whose maximal eigenvalue to estimate.\n\nMust be square and support premultiply (A*⋅).\n\niters::Int=1: Number of power iterations to run. (Recommended: iters ≤ 3)\n\nKeywords\n\np::Real=0.05: Probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains the maximal eigenvalue of A with probability 1 - p.\n\nReferences\n\nCorollary of Theorem 1 in [Dixon1983]\n\n\n\n"
},

{
    "location": "randomized.html#Extremal-eigenvalue-estimates-1",
    "page": "Randomized algorithms",
    "title": "Extremal eigenvalue estimates",
    "category": "section",
    "text": "reigmin\nreigmax"
},

{
    "location": "randomized.html#IterativeSolvers.rnorm",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.rnorm",
    "category": "Function",
    "text": "rnorm(A, mvps)\n\nCompute a probabilistic upper bound on the norm of a matrix A. ‖A‖ ≤ α √(2/π) maxᵢ ‖Aωᵢ‖ with probability p=α^(-mvps).\n\nArguments\n\nA: matrix whose norm to estimate.\nmvps::Int: number of matrix-vector products to compute.\n\nKeywords\n\np::Real=0.05: probability of upper bound failing.\n\nOutput\n\nEstimate of ‖A‖.\n\nSee also rnorms for a different estimator that uses  premultiplying by both A and A'.\n\nReferences\n\nLemma 4.1 of Halko2011\n\n\n\n"
},

{
    "location": "randomized.html#IterativeSolvers.rnorms",
    "page": "Randomized algorithms",
    "title": "IterativeSolvers.rnorms",
    "category": "Function",
    "text": "rnorms(A, iters=1)\n\nEstimate matrix norm randomly using A'A.\n\nCompute a probabilistic upper bound on the norm of a matrix A.\n\nρ = √(‖(A'A)ʲω‖/‖(A'A)ʲ⁻¹ω‖)\n\nwhich is an estimate of the spectral norm of A produced by iters steps of the power method starting with normalized ω, is a lower bound on the true norm by a factor\n\nρ ≤ α ‖A‖\n\nwith probability greater than 1 - p, where p = 4\\sqrt(n/(iters-1)) α^(-2iters).\n\nArguments\n\nA: matrix whose norm to estimate.\niters::Int = 1: mumber of power iterations to perform.\n\nKeywords\n\np::Real = 0.05: probability of upper bound failing.\nAt = A': Transpose of A.\n\nOutput\n\nEstimate of ‖A‖.\n\nSee also rnorm for a different estimator that does not require premultiplying by A'\n\nReferences\n\nAppendix of [Liberty2007].\n\n\n\n"
},

{
    "location": "randomized.html#Norm-estimate-1",
    "page": "Randomized algorithms",
    "title": "Norm estimate",
    "category": "section",
    "text": "rnorm\nrnorms[Halko2011]: Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. \"Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions.\" SIAM review 53.2 (2011): 217-288.[Dixon1983]: Dixon, John D. \"Estimating extremal eigenvalues and condition numbers of matrices.\" SIAM Journal on Numerical Analysis 20.4 (1983): 812-814.[Liberty2007]: Liberty, Edo, et al. \"Randomized algorithms for the low-rank approximation of matrices.\" Proceedings of the National Academy of Sciences 104.51 (2007): 20167-20172."
},

{
    "location": "iterators.html#",
    "page": "The iterator approach",
    "title": "The iterator approach",
    "category": "page",
    "text": ""
},

{
    "location": "iterators.html#Iterators-1",
    "page": "The iterator approach",
    "title": "Iterative solvers as iterators",
    "category": "section",
    "text": "In advanced use cases you might want to access the internal data structures of the solver, inject code to be run after each iteration, have total control over allocations or reduce overhead in initialization. The iterator approach of IterativeSolvers.jl makes this possible.note: Note\nAt this point BiCGStab(l), CG, Chebyshev, GMRES, MINRES and the stationary methods are implemented as iterators. However, the package does not yet export the iterators and helper methods themselves."
},

{
    "location": "iterators.html#How-iterators-are-implemented-1",
    "page": "The iterator approach",
    "title": "How iterators are implemented",
    "category": "section",
    "text": "The solvers listed above are basically a thin wrapper around an iterator. Among other things, they initialize the iterable, loop through the iterator and return the result:function my_solver!(x, A, b)\n    iterable = MySolverIterable(x, A, b)\n    for item in iterable end\n    return iterable.x\nendRather than calling my_solver!(x, A, b), you could also initialize the iterable yourself and perform the for loop."
},

{
    "location": "iterators.html#Example:-avoiding-unnecessary-initialization-1",
    "page": "The iterator approach",
    "title": "Example: avoiding unnecessary initialization",
    "category": "section",
    "text": "The Jacobi method for SparseMatrixCSC has some overhead in intialization; not only do we need to allocate a temporary vector, we also have to search for indices of the diagonal (and check their values are nonzero). The current implementation initializes the iterable as:jacobi_iterable(x, A::SparseMatrixCSC, b; maxiter::Int = 10) =\n    JacobiIterable{eltype(x), typeof(x)}(OffDiagonal(A, DiagonalIndices(A)), x, similar(x), b, maxiter)Now if you apply Jacobi iteration multiple times with the same matrix for just a few iterations, it makes sense to initialize the iterable only once and reuse it afterwards:A = sprand(10_000, 10_000, 10 / 10_000) + 20I\nb1 = rand(10_000)\nb2 = rand(10_000)\nx = rand(10_000)\n\nmy_iterable = IterativeSolvers.jacobi_iterable(x, A, b1, maxiter = 4)\n\nfor item in my_iterable \n    println(\"Iteration for rhs 1\")\nend\n\n@show norm(b1 - A * x) / norm(b1)\n\n# Copy the next right-hand side into the iterable\ncopy!(my_iterable.b, b2)\n\nfor item in my_iterable\n    println(\"Iteration for rhs 2\")\nend\n\n@show norm(b2 - A * x) / norm(b2)This would output:Iteration for rhs 1\nIteration for rhs 1\nIteration for rhs 1\nIteration for rhs 1\nnorm(b1 - A * x) / norm(b1) = 0.08388528015119746\nIteration for rhs 2\nIteration for rhs 2\nIteration for rhs 2\nIteration for rhs 2\nnorm(b2 - A * x) / norm(b2) = 0.0003681972775644809"
},

{
    "location": "iterators.html#Other-use-cases-1",
    "page": "The iterator approach",
    "title": "Other use cases",
    "category": "section",
    "text": "Other use cases include: computing the (harmonic) Ritz values from the Hessenberg matrix in GMRES;\ncomparing the approximate residual of methods such as GMRES and BiCGStab(l) with the true residual during the iterations;\nupdating a preconditioner in flexible methods."
},

{
    "location": "about/CONTRIBUTING.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "about/CONTRIBUTING.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "Contributions are always welcome, as are feature requests and suggestions. Feel free to open an issue or pull request at any time.It is important to note that almost every method in the package has documentation, to know what it does simply use ?<method> in the terminal.julia> using IterativeSolvers\n\nhelp?> IterativeSolvers.Adivtype\n  Adivtype(A, b)\n\n  Determine type of the division of an element of b against an element of A:\n\n  typeof(one(eltype(b))/one(eltype(A)))"
},

{
    "location": "about/CONTRIBUTING.html#Setting-workspace-up-1",
    "page": "Contributing",
    "title": "Setting workspace up",
    "category": "section",
    "text": "Julia's internal package manager makes it easy to install and modify packages from Github. Any package hosted on Github can be installed via Pkg.clone by providing the repository's URL, so installing a fork on your system is a simple task.Pkg.clone(\"https://github.com/johndoe/IterativeSolvers.jl\")It is to note here if you have the original package installed the fork will replace it, this is not a problem.Now find your fork's location.Pkg.dir(\"IterativeSolvers\")Once there you will notice you are on the master branch, whenever a package is imported Julia will use the code in the current branch, this means checking out other git branches will let you use/test whatever there is."
},

{
    "location": "about/CONTRIBUTING.html#Adding-or-modifying-iterative-methods-1",
    "page": "Contributing",
    "title": "Adding or modifying iterative methods",
    "category": "section",
    "text": "Each iterative method method must log information using the inner ConvergenceHistory type.There are two types of ConvergenceHistory: plain and restarted. The only real difference between the two is how they are plotted and how the number of restarts is calculated, everything else is the same.Before logging information space must always be reserved.log = ConvergenceHistory()\nlog[:tol] = tol\nreserve!(log,:betas, maxiter) # Vector of length maxiter\nreserve!(log,:conv, maxiter, T=BitArray) # Vector of length maxiter\nreserve!(log,:ritz, maxiter, k) # Matrix of size (maxiter, k)To store information at each iteration use push!.push!(log, :conv, conv)\npush!(log, :ritz, F[:S][1:k])\npush!(log, :betas, L.β)To advance the log index to the next iteration use nextiter!.nextiter!(log)A more detailed explanation of all the functions is in both the public and internal documentation of ConvergenceHistory.The most rich example of the usage of ConvergenceHistory is in svdl."
},

{
    "location": "about/license.html#",
    "page": "License",
    "title": "License",
    "category": "page",
    "text": ""
},

{
    "location": "about/license.html#License-(MIT)-1",
    "page": "License",
    "title": "License (MIT)",
    "category": "section",
    "text": "Copyright (c) 2013--2016 The Julia Language\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of\nthis software and associated documentation files (the \"Software\"), to deal in\nthe Software without restriction, including without limitation the rights to\nuse, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of\nthe Software, and to permit persons to whom the Software is furnished to do so,\nsubject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS\nFOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR\nCOPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER\nIN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN\nCONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

]}
