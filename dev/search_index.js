var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#IterativeSolvers.jl-1",
    "page": "Home",
    "title": "IterativeSolvers.jl",
    "category": "section",
    "text": "IterativeSolvers.jl is a Julia package that provides efficient iterative algorithms for solving large linear systems, eigenproblems, and singular value problems. Most of the methods can be used matrix-free.For bug reports, feature requests and questions please submit an issue. If you\'re interested in contributing, please see the Contributing guide.For more information on future methods have a look at the package roadmap."
},

{
    "location": "#What-method-should-I-use-for-linear-systems?-1",
    "page": "Home",
    "title": "What method should I use for linear systems?",
    "category": "section",
    "text": "When solving linear systems Ax = b for a square matrix A there are quite some options. The typical choices are listed below:Method When to use it\nConjugate Gradients Best choice for symmetric, positive-definite matrices\nMINRES For symmetric, indefinite matrices\nGMRES For nonsymmetric matrices when a good preconditioner is available\nIDR(s) For nonsymmetric, strongly indefinite problems without a good preconditioner\nBiCGStab(l) Otherwise for nonsymmetric problemsWe also offer Chebyshev iteration as an alternative to Conjugate Gradients when bounds on the spectrum are known.Stationary methods like Jacobi, Gauss-Seidel, SOR and SSOR can be used as smoothers to reduce high-frequency components in the error in just a few iterations.When solving least-squares problems we currently offer just LSMR and LSQR."
},

{
    "location": "#Eigenproblems-and-SVD-1",
    "page": "Home",
    "title": "Eigenproblems and SVD",
    "category": "section",
    "text": "For the Singular Value Decomposition we offer SVDL, which is the Golub-Kahan-Lanczos procedure.For eigenvalue problems we have at this point just the Power Method and some convenience wrappers to do shift-and-invert."
},

{
    "location": "getting_started/#",
    "page": "Getting started",
    "title": "Getting started",
    "category": "page",
    "text": ""
},

{
    "location": "getting_started/#Getting-started-1",
    "page": "Getting started",
    "title": "Getting started",
    "category": "section",
    "text": ""
},

{
    "location": "getting_started/#Installation-1",
    "page": "Getting started",
    "title": "Installation",
    "category": "section",
    "text": "The package can be installed via Julia\'s package manager.julia> Pkg.add(\"IterativeSolvers\")"
},

{
    "location": "getting_started/#Interface-1",
    "page": "Getting started",
    "title": "Interface",
    "category": "section",
    "text": "Virtually all solvers have the common function declarations:solver(A, args...; kwargs...)\nsolver!(x, A, args...; kwargs...)where A is a linear operator and x an initial guess. The second declaration also updates x in-place."
},

{
    "location": "getting_started/#matrixfree-1",
    "page": "Getting started",
    "title": "Explicit matrices and the matrix-free approach",
    "category": "section",
    "text": "Rather than constructing an explicit matrix A of the type Matrix or SparseMatrixCSC, it is also possible to pass a general linear operator that performs matrix operations implicitly. This is called the matrix-free approach.For matrix-free types of A the following interface is expected to be defined:A*v computes the matrix-vector product on a v::AbstractVector;\nmul!(y, A, v) computes the matrix-vector product on a v::AbstractVector in-place;\neltype(A) returns the element type implicit in the equivalent matrix representation of A;\nsize(A, d) returns the nominal dimensions along the dth axis in the equivalent matrix representation of A.tip: Matrix-free with LinearMaps.jl\nWe strongly recommend LinearMaps.jl for matrix-free linear operators, as it implements the above methods already for you; you just have to write the action of the linear map."
},

{
    "location": "getting_started/#Additional-arguments-1",
    "page": "Getting started",
    "title": "Additional arguments",
    "category": "section",
    "text": "Keyword names will vary depending on the method, however some of them will always have the same spelling:tol: (relative) stopping tolerance of the method;\nverbose: print information during the iterations;\nmaxiter: maximum number of allowed iterations;\nPl and Pr: left and right preconditioner. See Preconditioning;\nlog::Bool = false: output an extra element of type ConvergenceHistory containing the convergence history."
},

{
    "location": "getting_started/#log-keyword-1",
    "page": "Getting started",
    "title": "log keyword",
    "category": "section",
    "text": "Most solvers contain the log keyword. This is to be used when obtaining more information is required, to use it place the set log to true.x, ch = cg(Master, rand(10, 10), rand(10) log=true)\nsvd, L, ch = svdl(Master, rand(100, 100), log=true)The function will now return one more parameter of type ConvergenceHistory."
},

{
    "location": "getting_started/#IterativeSolvers.ConvergenceHistory",
    "page": "Getting started",
    "title": "IterativeSolvers.ConvergenceHistory",
    "category": "type",
    "text": "Store general and in-depth information about an iterative method.\n\nFields\n\nmvps::Int: number of matrix vector products.\n\nmtvps::Int: number of transposed matrix-vector products\n\niters::Int: iterations taken by the method.\n\nrestart::T: restart relevant information.\n\nT == Int: iterations per restart.\nT == Nothing: methods without restarts.\n\nisconverged::Bool: convergence of the method.\n\ndata::Dict{Symbol,Any}: Stores all the information stored during the method execution. It stores tolerances, residuals and other information, e.g. Ritz values in svdl.\n\nConstructors\n\nConvergenceHistory()\nConvergenceHistory(restart)\n\nCreate ConvergenceHistory with empty fields.\n\nArguments\n\nrestart: number of iterations per restart.\n\nPlots\n\nSupports plots using the Plots.jl package via a type recipe. Vectors are ploted as series and matrices as scatterplots.\n\nImplements\n\nBase: getindex, setindex!, push!\n\n\n\n\n\n"
},

{
    "location": "getting_started/#ConvergenceHistory-1",
    "page": "Getting started",
    "title": "ConvergenceHistory",
    "category": "section",
    "text": "A ConvergenceHistory instance stores information of a solver.Number of iterations.ch.itersConvergence status.ch.isconvergedStopping tolerances. (A Symbol key is needed to access)ch[:tol]Maximum number of iterations per restart. (Only on restarted methods)nrests(ch)Number of matrix-vectors and matrix-transposed-vector products.nprods(ch)Data stored on each iteration, accessed information can be either a vector or matrix. This data can be a lot of things, most commonly residual. (A Symbol key is needed to access)ch[:resnorm] #Vector or Matrix\nch[:resnorm, x] #Vector or Matrix element\nch[:resnorm, x, y] #Matrix elementConvergenceHistory"
},

{
    "location": "getting_started/#Plotting-1",
    "page": "Getting started",
    "title": "Plotting",
    "category": "section",
    "text": "ConvergeHistory provides a recipe to use with the package Plots.jl, this makes it really easy to plot on different plot backends. There are two recipes provided:One for the whole ConvergenceHistory.plot(ch)The other one to plot data binded to a key._, ch = gmres(rand(10,10), rand(10), maxiter = 100, log=true)\nplot(ch, :resnorm, sep = :blue)Plot additional keywordssep::Symbol = :white: color of the line separator in restarted methods."
},

{
    "location": "preconditioning/#",
    "page": "Preconditioning",
    "title": "Preconditioning",
    "category": "page",
    "text": ""
},

{
    "location": "preconditioning/#Preconditioning-1",
    "page": "Preconditioning",
    "title": "Preconditioning",
    "category": "section",
    "text": "Many iterative solvers have the option to provide left and right preconditioners (Pl and Pr resp.) in order to speed up convergence or prevent stagnation. They transform a problem Ax = b into a better conditioned system (P_l^-1AP_r^-1)y = P_l^-1b, where x = P_r^-1y.These preconditioners should support the operationsldiv!(y, P, x) computes P \\ x in-place of y;\nldiv!(P, x) computes P \\ x in-place of x;\nand P \\ x.If no preconditioners are passed to the solver, the method will default toPl = Pr = IterativeSolvers.Identity()"
},

{
    "location": "preconditioning/#Available-preconditioners-1",
    "page": "Preconditioning",
    "title": "Available preconditioners",
    "category": "section",
    "text": "IterativeSolvers.jl itself does not provide any other preconditioners besides Identity(), but recommends the following external packages:ILU.jl for incomplete LU decompositions (using drop tolerance);\nIncompleteSelectedInversion.jl for incomplete LDLt decompositions.\nAMG.jl for some algebraic multigrid (AMG) preconditioners.\nPreconditioners.jl which wraps a bunch of preconditioners from other packages. If you are a beginner or want to try different ones quickly, this is good starting place."
},

{
    "location": "linear_systems/cg/#",
    "page": "Conjugate Gradients",
    "title": "Conjugate Gradients",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/cg/#CG-1",
    "page": "Conjugate Gradients",
    "title": "Conjugate Gradients (CG)",
    "category": "section",
    "text": "Conjugate Gradients solves Ax = b approximately for x where A is a symmetric, positive-definite linear operator and b the right-hand side vector. The method uses short recurrences and therefore has fixed memory costs and fixed computational costs per iteration."
},

{
    "location": "linear_systems/cg/#IterativeSolvers.cg",
    "page": "Conjugate Gradients",
    "title": "IterativeSolvers.cg",
    "category": "function",
    "text": "cg(A, b; kwargs...) -> x, [history]\n\nSame as cg!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/cg/#IterativeSolvers.cg!",
    "page": "Conjugate Gradients",
    "title": "IterativeSolvers.cg!",
    "category": "function",
    "text": "cg!(x, A, b; kwargs...) -> x, [history]\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\nstatevars::CGStateVariables: Has 3 arrays similar to x to hold intermediate results;\ninitially_zero::Bool: If true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;\nPl = Identity(): left preconditioner of the method. Should be symmetric, positive-definite like A;\ntol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol;\nmaxiter::Int = size(A,2): maximum number of iterations;\nverbose::Bool = false: print method information;\nlog::Bool = false: keep track of the residual norm in each iteration.\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/cg/#Usage-1",
    "page": "Conjugate Gradients",
    "title": "Usage",
    "category": "section",
    "text": "cg\ncg!"
},

{
    "location": "linear_systems/cg/#On-the-GPU-1",
    "page": "Conjugate Gradients",
    "title": "On the GPU",
    "category": "section",
    "text": "The method should work fine on the GPU. As a minimal working example, consider:using LinearAlgebra, CuArrays, IterativeSolvers\n\nn = 100\nA = cu(rand(n, n))\nA = A + A\' + 2*n*I\nb = cu(rand(n))\nx = cg(A, b)note: Note\nMake sure that all state vectors are stored on the GPU. For instance when calling cg!(x, A, b), one might have an issue when x is stored on the GPU, while b is stored on the CPU – IterativeSolvers.jl does not copy the vectors to the same device."
},

{
    "location": "linear_systems/cg/#Implementation-details-1",
    "page": "Conjugate Gradients",
    "title": "Implementation details",
    "category": "section",
    "text": "The current implementation follows a rather standard approach. Note that preconditioned CG (or PCG) is slightly different from ordinary CG, because the former must compute the residual explicitly, while it is available as byproduct in the latter. Our implementation of CG ensures the minimal number of vector operations.tip: Tip\nCG can be used as an iterator."
},

{
    "location": "linear_systems/chebyshev/#",
    "page": "Chebyshev iteration",
    "title": "Chebyshev iteration",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/chebyshev/#Chebyshev-1",
    "page": "Chebyshev iteration",
    "title": "Chebyshev iteration",
    "category": "section",
    "text": "Chebyshev iteration solves the problem Ax=b approximately for x where A is a symmetric, definite linear operator and b the right-hand side vector. The methods assumes the interval lambda_min lambda_max containing all eigenvalues of A is known, so that x can be iteratively constructed via a Chebyshev polynomial with zeros in this interval. This polynomial ultimately acts as a filter that removes components in the direction of the eigenvectors from the initial residual.The main advantage with respect to Conjugate Gradients is that BLAS1 operations such as inner products are avoided."
},

{
    "location": "linear_systems/chebyshev/#IterativeSolvers.chebyshev",
    "page": "Chebyshev iteration",
    "title": "IterativeSolvers.chebyshev",
    "category": "function",
    "text": "chebyshev(A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]\n\nSame as chebyshev!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/chebyshev/#IterativeSolvers.chebyshev!",
    "page": "Chebyshev iteration",
    "title": "IterativeSolvers.chebyshev!",
    "category": "function",
    "text": "chebyshev!(x, A, b, λmin::Real, λmax::Real; kwargs...) -> x, [history]\n\nSolve Ax = b for symmetric, definite matrices A using Chebyshev iteration.\n\nArguments\n\nx: initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side;\nλmin::Real: lower bound for the real eigenvalues\nλmax::Real: upper bound for the real eigenvalues\n\nKeywords\n\ninitially_zero::Bool = false: if true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;\ntol: tolerance for stopping condition |r_k| / |r_0| ≤ tol.\nmaxiter::Int = size(A, 2): maximum number of inner iterations of GMRES;\nPl = Identity(): left preconditioner;\nlog::Bool = false: keep track of the residual norm in each iteration;\nverbose::Bool = false: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/chebyshev/#Usage-1",
    "page": "Chebyshev iteration",
    "title": "Usage",
    "category": "section",
    "text": "chebyshev\nchebyshev!"
},

{
    "location": "linear_systems/chebyshev/#Implementation-details-1",
    "page": "Chebyshev iteration",
    "title": "Implementation details",
    "category": "section",
    "text": "warning: BLAS1 operations\nAlthough the method is often used to avoid computation of inner products, the stopping criterion is still based on the residual norm. Hence the current implementation is not free of BLAS1 operations.tip: Tip\nChebyshev iteration can be used as an iterator."
},

{
    "location": "linear_systems/minres/#",
    "page": "MINRES",
    "title": "MINRES",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/minres/#MINRES-1",
    "page": "MINRES",
    "title": "MINRES",
    "category": "section",
    "text": "MINRES is a short-recurrence version of GMRES for solving Ax = b approximately for x where A is a symmetric, Hermitian, skew-symmetric or skew-Hermitian linear operator and b the right-hand side vector."
},

{
    "location": "linear_systems/minres/#IterativeSolvers.minres",
    "page": "MINRES",
    "title": "IterativeSolvers.minres",
    "category": "function",
    "text": "minres(A, b; kwargs...) -> x, [history]\n\nSame as minres!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/minres/#IterativeSolvers.minres!",
    "page": "MINRES",
    "title": "IterativeSolvers.minres!",
    "category": "function",
    "text": "minres!(x, A, b; kwargs...) -> x, [history]\n\nSolve Ax = b for (skew-)Hermitian matrices A using MINRES.\n\nArguments\n\nx: initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool = false: if true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;\nskew_hermitian::Bool = false: if true assumes that A is skew-symmetric or skew-Hermitian;\ntol: tolerance for stopping condition |r_k| / |r_0| ≤ tol. Note that the residual is computed only approximately;\nmaxiter::Int = size(A, 2): maximum number of iterations;\nPl: left preconditioner;\nPr: right preconditioner;\nlog::Bool = false: keep track of the residual norm in each iteration;\nverbose::Bool = false: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/minres/#Usage-1",
    "page": "MINRES",
    "title": "Usage",
    "category": "section",
    "text": "minres\nminres!"
},

{
    "location": "linear_systems/minres/#Implementation-details-1",
    "page": "MINRES",
    "title": "Implementation details",
    "category": "section",
    "text": "MINRES exploits the tridiagonal structure of the Hessenberg matrix. Although MINRES is mathematically equivalent to GMRES, it might not be equivalent in finite precision. MINRES updates the solution asx = x_0 + (V R^-1) (Q^*r_0e_1)where V is the orthonormal basis for the Krylov subspace and QR is the QR-decomposition of the Hessenberg matrix. Note that the brackets are placed slightly differently from how GMRES would update the residual.MINRES computes V and W = VR^-1 via a three-term recurrence, using only the last column of R Therefore we pre-allocate only six vectors, save only the last two entries of Q^*r_0e_1 and part of the last column of the Hessenberg matrix.note: Real and complex arithmetic\nIf A is Hermitian, then the Hessenberg matrix will be real. This is exploited in the current implementation.If A is skew-Hermitian, the diagonal of the Hessenberg matrix will be imaginary, and hence we use complex arithmetic in that case.tip: Tip\nMINRES can be used as an iterator."
},

{
    "location": "linear_systems/bicgstabl/#",
    "page": "BiCGStab(l)",
    "title": "BiCGStab(l)",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/bicgstabl/#BiCGStabl-1",
    "page": "BiCGStab(l)",
    "title": "BiCGStab(l)",
    "category": "section",
    "text": "BiCGStab(l) solves the problem Ax = b approximately for x where A is a general, linear operator and b the right-hand side vector. The methods combines BiCG with l GMRES iterations, resulting in a short-reccurence iteration. As a result the memory is fixed as well as the computational costs per iteration."
},

{
    "location": "linear_systems/bicgstabl/#IterativeSolvers.bicgstabl",
    "page": "BiCGStab(l)",
    "title": "IterativeSolvers.bicgstabl",
    "category": "function",
    "text": "bicgstabl(A, b, l; kwargs...) -> x, [history]\n\nSame as bicgstabl!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/bicgstabl/#IterativeSolvers.bicgstabl!",
    "page": "BiCGStab(l)",
    "title": "IterativeSolvers.bicgstabl!",
    "category": "function",
    "text": "bicgstabl!(x, A, b, l; kwargs...) -> x, [history]\n\nArguments\n\nA: linear operator;\nb: right hand side (vector);\nl::Int = 2: Number of GMRES steps.\n\nKeywords\n\nmax_mv_products::Int = size(A, 2): maximum number of matrix vector products.\n\nFor BiCGStab(l) this is a less dubious term than \"number of iterations\";\n\nPl = Identity(): left preconditioner of the method;\ntol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol.  Note that (1) the true residual norm is never computed during the iterations,  only an approximation; and (2) if a preconditioner is given, the stopping condition is based on the  preconditioned residual.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/bicgstabl/#Usage-1",
    "page": "BiCGStab(l)",
    "title": "Usage",
    "category": "section",
    "text": "bicgstabl\nbicgstabl!"
},

{
    "location": "linear_systems/bicgstabl/#Implementation-details-1",
    "page": "BiCGStab(l)",
    "title": "Implementation details",
    "category": "section",
    "text": "The method is based on the original article [Sleijpen1993], but does not implement later improvements. The normal equations arising from the GMRES steps are solved without orthogonalization. Hence the method should only be reliable for relatively small values of l.The r and u factors are pre-allocated as matrices of size n times (l + 1), so that BLAS2 methods can be used. Also the random shadow residual is pre-allocated as a vector. Hence the storage costs are approximately 2l + 3 vectors.tip: Tip\nBiCGStabl(l) can be used as an iterator.[Sleijpen1993]: Sleijpen, Gerard LG, and Diederik R. Fokkema. \"BiCGstab(l) for  linear equations involving unsymmetric matrices with complex spectrum.\"  Electronic Transactions on Numerical Analysis 1.11 (1993): 2000."
},

{
    "location": "linear_systems/idrs/#",
    "page": "IDR(s)",
    "title": "IDR(s)",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/idrs/#IDRs-1",
    "page": "IDR(s)",
    "title": "IDR(s)",
    "category": "section",
    "text": "The Induced Dimension Reduction method is a family of simple and fast Krylov subspace algorithms for solving large nonsymmetric linear systems. The idea behind the IDR(s) variant is to generate residuals that are in the nested subspaces of shrinking dimensions."
},

{
    "location": "linear_systems/idrs/#IterativeSolvers.idrs",
    "page": "IDR(s)",
    "title": "IterativeSolvers.idrs",
    "category": "function",
    "text": "idrs(A, b; s = 8) -> x, [history]\n\nSame as idrs!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/idrs/#IterativeSolvers.idrs!",
    "page": "IDR(s)",
    "title": "IterativeSolvers.idrs!",
    "category": "function",
    "text": "idrs!(x, A, b; s = 8) -> x, [history]\n\nSolve the problem Ax = b approximately with IDR(s), where s is the dimension of the shadow space.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ns::Integer = 8: dimension of the shadow space;\ntol: relative tolerance;\nmaxiter::Int = size(A, 2): maximum number of iterations;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/idrs/#Usage-1",
    "page": "IDR(s)",
    "title": "Usage",
    "category": "section",
    "text": "idrs\nidrs!"
},

{
    "location": "linear_systems/idrs/#Implementation-details-1",
    "page": "IDR(s)",
    "title": "Implementation details",
    "category": "section",
    "text": "The current implementation is based on the MATLAB version by Van Gijzen and Sonneveld. For background see [Sonneveld2008], [VanGijzen2011] and the IDR(s) webpage.[Sonneveld2008]: IDR(s): a family of simple and fast algorithms for solving large nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035–1062, 2008[VanGijzen2011]: Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011"
},

{
    "location": "linear_systems/gmres/#",
    "page": "Restarted GMRES",
    "title": "Restarted GMRES",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/gmres/#GMRES-1",
    "page": "Restarted GMRES",
    "title": "Restarted GMRES",
    "category": "section",
    "text": "GMRES solves the problem Ax = b approximately for x where A is a general, linear operator and b the right-hand side vector. The method is optimal in the sense that it selects the solution with minimal residual from a Krylov subspace, but the price of optimality is increasing storage and computation effort per iteration. Restarts are necessary to fix these costs."
},

{
    "location": "linear_systems/gmres/#IterativeSolvers.gmres",
    "page": "Restarted GMRES",
    "title": "IterativeSolvers.gmres",
    "category": "function",
    "text": "gmres(A, b; kwargs...) -> x, [history]\n\nSame as gmres!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/gmres/#IterativeSolvers.gmres!",
    "page": "Restarted GMRES",
    "title": "IterativeSolvers.gmres!",
    "category": "function",
    "text": "gmres!(x, A, b; kwargs...) -> x, [history]\n\nSolves the problem Ax = b with restarted GMRES.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ninitially_zero::Bool: If true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;\ntol: relative tolerance;\nrestart::Int = min(20, size(A, 2)): restarts GMRES after specified number of iterations;\nmaxiter::Int = size(A, 2): maximum number of inner iterations of GMRES;\nPl: left preconditioner;\nPr: right preconditioner;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximate solution.\n\nif log is true\n\nx: approximate solution;\nhistory: convergence history.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/gmres/#Usage-1",
    "page": "Restarted GMRES",
    "title": "Usage",
    "category": "section",
    "text": "gmres\ngmres!"
},

{
    "location": "linear_systems/gmres/#Implementation-details-1",
    "page": "Restarted GMRES",
    "title": "Implementation details",
    "category": "section",
    "text": "The implementation pre-allocates a matrix V of size n by restart whose columns form an orthonormal basis for the Krylov subspace. This allows BLAS2 operations when updating the solution vector x. The Hessenberg matrix is also pre-allocated.Modified Gram-Schmidt is used to orthogonalize the columns of V.The computation of the residual norm is implemented in a non-standard way, namely keeping track of a vector gamma in the null-space of H_k^*, which is the adjoint of the (k + 1) times k Hessenberg matrix H_k at the kth iteration. Only when x needs to be updated is the Hessenberg matrix mutated with Givens rotations.tip: Tip\nGMRES can be used as an iterator. This makes it possible to access the Hessenberg matrix and Krylov basis vectors during the iterations."
},

{
    "location": "linear_systems/lsmr/#",
    "page": "LSMR",
    "title": "LSMR",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/lsmr/#LSMR-1",
    "page": "LSMR",
    "title": "LSMR",
    "category": "section",
    "text": "Least-squares minimal residual"
},

{
    "location": "linear_systems/lsmr/#IterativeSolvers.lsmr",
    "page": "LSMR",
    "title": "IterativeSolvers.lsmr",
    "category": "function",
    "text": "lsmr(A, b; kwrags...) -> x, [history]\n\nSame as lsmr!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/lsmr/#IterativeSolvers.lsmr!",
    "page": "LSMR",
    "title": "IterativeSolvers.lsmr!",
    "category": "function",
    "text": "lsmr!(x, A, b; kwargs...) -> x, [history]\n\nMinimizes Ax - b^2 + λx^2 in the Euclidean norm. If multiple solutions exists the minimum norm solution is returned.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying MINRES to the normal equations (A^*A + λ^2I)x = A^*b, but has better numerical properties, especially if A is ill-conditioned.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\nconlim::Number = 1e8: stopping tolerance. lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting\natol = btol = conlim = zero, but the number of iterations may then be excessive.\nmaxiter::Int = maximum(size(A)): maximum number of iterations.\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nReturn values\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n:btol => ::Real: btol stopping tolerance.\n:ctol => ::Real: ctol stopping tolerance.\n:anorm => ::Real: anorm.\n:rnorm => ::Real: rnorm.\n:cnorm => ::Real: cnorm.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/lsmr/#Usage-1",
    "page": "LSMR",
    "title": "Usage",
    "category": "section",
    "text": "lsmr\nlsmr!"
},

{
    "location": "linear_systems/lsmr/#Implementation-details-1",
    "page": "LSMR",
    "title": "Implementation details",
    "category": "section",
    "text": "Adapted from: http://web.stanford.edu/group/SOL/software/lsmr/"
},

{
    "location": "linear_systems/lsqr/#",
    "page": "LSQR",
    "title": "LSQR",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/lsqr/#LSQR-1",
    "page": "LSQR",
    "title": "LSQR",
    "category": "section",
    "text": ""
},

{
    "location": "linear_systems/lsqr/#IterativeSolvers.lsqr",
    "page": "LSQR",
    "title": "IterativeSolvers.lsqr",
    "category": "function",
    "text": "lsqr(A, b; kwrags...) -> x, [history]\n\nSame as lsqr!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/lsqr/#IterativeSolvers.lsqr!",
    "page": "LSQR",
    "title": "IterativeSolvers.lsqr!",
    "category": "function",
    "text": "lsqr!(x, A, b; kwargs...) -> x, [history]\n\nMinimizes Ax - b^2 + damp*x^2 in the Euclidean norm. If multiple solutions exists returns the minimal norm solution.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying CG to the normal equations (A^*A + λ^2I)x = A^*b but has better numerical properties, especially if A is ill-conditioned.\n\nArguments\n\nx: Initial guess, will be updated in-place;\nA: linear operator;\nb: right-hand side.\n\nKeywords\n\ndamp::Number = 0: damping parameter.\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\nmaxiter::Int = maximum(size(A)): maximum number of iterations.\nverbose::Bool = false: print method information.\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nReturn values\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n:btol => ::Real: btol stopping tolerance.\n:ctol => ::Real: ctol stopping tolerance.\n:anorm => ::Real: anorm.\n:rnorm => ::Real: rnorm.\n:cnorm => ::Real: cnorm.\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/lsqr/#Usage-1",
    "page": "LSQR",
    "title": "Usage",
    "category": "section",
    "text": "lsqr\nlsqr!"
},

{
    "location": "linear_systems/lsqr/#Implementation-details-1",
    "page": "LSQR",
    "title": "Implementation details",
    "category": "section",
    "text": "Adapted from: http://web.stanford.edu/group/SOL/software/lsqr/."
},

{
    "location": "linear_systems/stationary/#",
    "page": "Stationary methods",
    "title": "Stationary methods",
    "category": "page",
    "text": ""
},

{
    "location": "linear_systems/stationary/#Stationary-1",
    "page": "Stationary methods",
    "title": "Stationary methods",
    "category": "section",
    "text": "Stationary methods are typically used as smoothers in multigrid methods, where only very few iterations are applied to get rid of high-frequency components in the error. The implementations of stationary methods have this goal in mind, which means there is no other stopping criterion besides the maximum number of iterations.note: CSC versus CSR\nJulia stores matrices column-major. In order to avoid cache misses, the implementations of our stationary methods traverse the matrices column-major. This deviates from classical textbook implementations. Also the SOR and SSOR methods cannot be computed efficiently in-place, but require a temporary vector.When it comes to SparseMatrixCSC, we precompute in all stationary methods an integer array of the indices of the diagonal to avoid expensive searches in each iteration."
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.jacobi",
    "page": "Stationary methods",
    "title": "IterativeSolvers.jacobi",
    "category": "function",
    "text": "jacobi(A, b) -> x\n\nSame as jacobi!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.jacobi!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.jacobi!",
    "category": "function",
    "text": "jacobi!(x, A::AbstractMatrix, b; maxiter=10) -> x\n\nPerforms exactly maxiter Jacobi iterations.\n\nAllocates a single temporary vector and traverses A columnwise.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\njacobi!(x, A::SparseMatrixCSC, b; maxiter=10) -> x\n\nPerforms exactly maxiter Jacobi iterations.\n\nAllocates a temporary vector and precomputes the diagonal indices.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#Jacobi-1",
    "page": "Stationary methods",
    "title": "Jacobi",
    "category": "section",
    "text": "jacobi\njacobi!"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.gauss_seidel",
    "page": "Stationary methods",
    "title": "IterativeSolvers.gauss_seidel",
    "category": "function",
    "text": "gauss_seidel(A, b) -> x\n\nSame as gauss_seidel!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.gauss_seidel!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.gauss_seidel!",
    "category": "function",
    "text": "gauss_seidel!(x, A::AbstractMatrix, b; maxiter=10) -> x\n\nPerforms exactly maxiter Gauss-Seidel iterations.\n\nWorks fully in-place and traverses A columnwise.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\ngauss_seidel!(x, A::SparseMatrixCSC, b; maxiter=10) -> x\n\nPerforms exactly maxiter Gauss-Seidel iterations.\n\nWorks fully in-place, but precomputes the diagonal indices.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#Gauss-Seidel-1",
    "page": "Stationary methods",
    "title": "Gauss-Seidel",
    "category": "section",
    "text": "gauss_seidel\ngauss_seidel!"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.sor",
    "page": "Stationary methods",
    "title": "IterativeSolvers.sor",
    "category": "function",
    "text": "sor(A, b, ω::Real) -> x\n\nSame as sor!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.sor!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.sor!",
    "category": "function",
    "text": "sor!(x, A::AbstractMatrix, b, ω::Real; maxiter=10) -> x\n\nPerforms exactly maxiter SOR iterations with relaxation parameter ω.\n\nAllocates a single temporary vector and traverses A columnwise.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\nsor!(x, A::SparseMatrixCSC, b, ω::Real; maxiter=10)\n\nPerforms exactly maxiter SOR iterations with relaxation parameter ω.\n\nAllocates a temporary vector and precomputes the diagonal indices.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#SOR-1",
    "page": "Stationary methods",
    "title": "Successive over-relaxation (SOR)",
    "category": "section",
    "text": "sor\nsor!"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.ssor",
    "page": "Stationary methods",
    "title": "IterativeSolvers.ssor",
    "category": "function",
    "text": "ssor(A, b, ω::Real) -> x\n\nSame as ssor!, but allocates a solution vector x initialized with zeros.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#IterativeSolvers.ssor!",
    "page": "Stationary methods",
    "title": "IterativeSolvers.ssor!",
    "category": "function",
    "text": "ssor!(x, A::AbstractMatrix, b, ω::Real; maxiter=10) -> x\n\nPerforms exactly maxiter SSOR iterations with relaxation parameter ω. Each iteration is basically a forward and backward sweep of SOR.\n\nAllocates a single temporary vector and traverses A columnwise.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\nssor!(x, A::SparseMatrixCSC, b, ω::Real; maxiter=10)\n\nPerforms exactly maxiter SSOR iterations with relaxation parameter ω. Each iteration is basically a forward and backward sweep of SOR.\n\nAllocates a temporary vector and precomputes the diagonal indices.\n\nThrows LinearAlgebra.SingularException when the diagonal has a zero. This check is performed once beforehand.\n\n\n\n\n\n"
},

{
    "location": "linear_systems/stationary/#SSOR-1",
    "page": "Stationary methods",
    "title": "Symmetric successive over-relaxation (SSOR)",
    "category": "section",
    "text": "ssor\nssor!tip: Tip\nAll stationary methods can be used a iterators."
},

{
    "location": "eigenproblems/power_method/#",
    "page": "Power method",
    "title": "Power method",
    "category": "page",
    "text": ""
},

{
    "location": "eigenproblems/power_method/#PowerMethod-1",
    "page": "Power method",
    "title": "(Inverse) power method",
    "category": "section",
    "text": "Solves the eigenproblem Ax = λx approximately where A is a general linear map. By default converges towards the dominant eigenpair (λ x) such that λ is largest. Shift-and-invert can be applied to target a specific eigenvalue near shift in the complex plane."
},

{
    "location": "eigenproblems/power_method/#IterativeSolvers.powm",
    "page": "Power method",
    "title": "IterativeSolvers.powm",
    "category": "function",
    "text": "powm(B; kwargs...) -> λ, x, [history]\n\nSee powm!. Calls powm!(B, x0; kwargs...) with x0 initialized as a random, complex unit vector.\n\n\n\n\n\n"
},

{
    "location": "eigenproblems/power_method/#IterativeSolvers.powm!",
    "page": "Power method",
    "title": "IterativeSolvers.powm!",
    "category": "function",
    "text": "powm!(B, x; shift = zero(eltype(B)), inverse::Bool = false, kwargs...) -> λ, x, [history]\n\nBy default finds the approximate eigenpair (λ, x) of B where |λ| is largest.\n\nArguments\n\nB: linear map, see the note below.\nx: normalized initial guess. Don\'t forget to use complex arithmetic when necessary.\n\nKeywords\n\ntol::Real = eps(real(eltype(B))) * size(B, 2) ^ 3: stopping tolerance for the residual norm;\nmaxiter::Integer = size(B,2): maximum number of iterations;\nlog::Bool: keep track of the residual norm in each iteration;\nverbose::Bool: print convergence information during the iterations.\n\nnote: Shift-and-invert\nWhen applying shift-and-invert to Ax = λx with invert = true and shift = ..., note that the role of B * b becomes computing inv(A - shift I) * b. So rather than passing the linear map A itself, pass a linear map B that has the action of shift-and-invert. The eigenvalue is transformed back to an eigenvalue of the actual matrix A.\n\nReturn values\n\nif log is false\n\nλ::Number approximate eigenvalue computed as the Rayleigh quotient;\nx::Vector approximate eigenvector.\n\nif log is true\n\nλ::Number: approximate eigenvalue computed as the Rayleigh quotient;\nx::Vector: approximate eigenvector;\nhistory: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance;\n:resnom => ::Vector: residual norm at each iteration.\n\nExamples\n\nusing LinearMaps\nσ = 1.0 + 1.3im\nA = rand(ComplexF64, 50, 50)\nF = lu(A - σ * I)\nFmap = LinearMap{ComplexF64}((y, x) -> ldiv!(y, F, x), 50, ismutating = true)\nλ, x = powm(Fmap, inverse = true, shift = σ, tol = 1e-4, maxiter = 200)\n\n\n\n\n\n"
},

{
    "location": "eigenproblems/power_method/#IterativeSolvers.invpowm",
    "page": "Power method",
    "title": "IterativeSolvers.invpowm",
    "category": "function",
    "text": "invpowm(B; shift = σ, kwargs...) -> λ, x, [history]\n\nFind the approximate eigenpair (λ, x) of A near shift, where B is a linear map that has the effect B * v = inv(A - σI) * v.\n\nThe method calls powm!(B, x0; inverse = true, shift = σ) with x0 a random, complex unit vector. See powm!\n\nExamples\n\nusing LinearMaps\nσ = 1.0 + 1.3im\nA = rand(ComplexF64, 50, 50)\nF = lu(A - σ * I)\nFmap = LinearMap{ComplexF64}((y, x) -> ldiv!(y, F, x), 50, ismutating = true)\nλ, x = invpowm(Fmap, shift = σ, tol = 1e-4, maxiter = 200)\n\n\n\n\n\n"
},

{
    "location": "eigenproblems/power_method/#IterativeSolvers.invpowm!",
    "page": "Power method",
    "title": "IterativeSolvers.invpowm!",
    "category": "function",
    "text": "invpowm!(B, x0; shift = σ, kwargs...) -> λ, x, [history]\n\nFind the approximate eigenpair (λ, x) of A near shift, where B is a linear map that has the effect B * v = inv(A - σI) * v.\n\nThe method calls powm!(B, x0; inverse = true, shift = σ). See powm!.\n\n\n\n\n\n"
},

{
    "location": "eigenproblems/power_method/#Usage-1",
    "page": "Power method",
    "title": "Usage",
    "category": "section",
    "text": "powm\npowm!\ninvpowm\ninvpowm!"
},

{
    "location": "eigenproblems/power_method/#Implementation-details-1",
    "page": "Power method",
    "title": "Implementation details",
    "category": "section",
    "text": "Storage requirements are 3 vectors: the approximate eigenvector x, the residual vector r and a temporary. The residual norm lags behind one iteration, as it is computed when Ax is performed. Therefore the final resdiual norm is even smaller."
},

{
    "location": "eigenproblems/lobpcg/#",
    "page": "LOBPCG",
    "title": "LOBPCG",
    "category": "page",
    "text": ""
},

{
    "location": "eigenproblems/lobpcg/#LOBPCG-1",
    "page": "LOBPCG",
    "title": "Locally optimal block preconditioned conjugate gradient (LOBPCG)",
    "category": "section",
    "text": "Solves the generalized eigenproblem Ax = λBx approximately where A and B are Hermitian linear maps, and B is positive definite. B is taken to be the identity by default. It can find the smallest (or largest) k eigenvalues and their corresponding eigenvectors which are B-orthonormal. It also admits a preconditioner and a \"constraints\" matrix C, such that the algorithm returns the smallest (or largest) eigenvalues associated with the eigenvectors in the nullspace of C\'B."
},

{
    "location": "eigenproblems/lobpcg/#IterativeSolvers.lobpcg",
    "page": "LOBPCG",
    "title": "IterativeSolvers.lobpcg",
    "category": "function",
    "text": "The Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG)\n\nFinds the nev extremal eigenvalues and their corresponding eigenvectors satisfying AX = λBX.\n\nA and B may be generic types but Base.mul!(C, AorB, X) must be defined for vectors and strided matrices X and C. size(A, i::Int) and eltype(A) must also be defined for A.\n\nlobpcg(A, [B,] largest, nev; kwargs...) -> results\n\nArguments\n\nA: linear operator;\nB: linear operator;\nlargest: true if largest eigenvalues are desired and false if smallest;\nnev: number of eigenvalues desired.\n\nKeywords\n\nlog::Bool: default is false; if true, results.trace will store iterations   states; if false only results.trace will be empty;\nP: preconditioner of residual vectors, must overload ldiv!;\nC: constraint to deflate the residual and solution vectors orthogonal   to a subspace; must overload mul!;\nmaxiter: maximum number of iterations; default is 200;\ntol::Real: tolerance to which residual vector norms must be under.\n\nOutput\n\nresults: a LOBPCGResults struct. r.λ and r.X store the eigenvalues and eigenvectors.\n\n\n\n\n\nlobpcg(A, [B,] largest, X0; kwargs...) -> results\n\nArguments\n\nA: linear operator;\nB: linear operator;\nlargest: true if largest eigenvalues are desired and false if smallest;\nX0: Initial guess, will not be modified. The number of columns is the number of eigenvectors desired.\n\nKeywords\n\nnot_zeros: default is false. If true, X0 will be assumed to not have any all-zeros column.\nlog::Bool: default is false; if true, results.trace will store iterations   states; if false only results.trace will be empty;\nP: preconditioner of residual vectors, must overload ldiv!;\nC: constraint to deflate the residual and solution vectors orthogonal   to a subspace; must overload mul!;\nmaxiter: maximum number of iterations; default is 200;\ntol::Real: tolerance to which residual vector norms must be under.\n\nOutput\n\nresults: a LOBPCGResults struct. r.λ and r.X store the eigenvalues and eigenvectors.\n\n\n\n\n\nlobpcg(A, [B,] largest, X0, nev; kwargs...) -> results\n\nArguments\n\nA: linear operator;\nB: linear operator;\nlargest: true if largest eigenvalues are desired and false if smallest;\nX0: block vectors such that the eigenvalues will be found size(X0, 2) at a time;   the columns are also used to initialize the first batch of Ritz vectors;\nnev: number of eigenvalues desired.\n\nKeywords\n\nlog::Bool: default is false; if true, results.trace will store iterations   states; if false only results.trace will be empty;\nP: preconditioner of residual vectors, must overload ldiv!;\nC: constraint to deflate the residual and solution vectors orthogonal   to a subspace; must overload mul!;\nmaxiter: maximum number of iterations; default is 200;\ntol::Real: tolerance to which residual vector norms must be under.\n\nOutput\n\nresults: a LOBPCGResults struct. r.λ and r.X store the eigenvalues and eigenvectors.\n\n\n\n\n\n"
},

{
    "location": "eigenproblems/lobpcg/#IterativeSolvers.lobpcg!",
    "page": "LOBPCG",
    "title": "IterativeSolvers.lobpcg!",
    "category": "function",
    "text": "lobpcg!(iterator::LOBPCGIterator; kwargs...) -> results\n\nArguments\n\niterator::LOBPCGIterator: a struct having all the variables required   for the LOBPCG algorithm.\n\nKeywords\n\nnot_zeros: default is false. If true, the initial Ritz vectors will be assumed to not have any all-zeros column.\nlog::Bool: default is false; if true, results.trace will store iterations   states; if false only results.trace will be empty;\nmaxiter: maximum number of iterations; default is 200;\ntol::Real: tolerance to which residual vector norms must be under.\n\nOutput\n\nresults: a LOBPCGResults struct. r.λ and r.X store the eigenvalues and eigenvectors.\n\n\n\n\n\n"
},

{
    "location": "eigenproblems/lobpcg/#Usage-1",
    "page": "LOBPCG",
    "title": "Usage",
    "category": "section",
    "text": "lobpcg\nlobpcg!"
},

{
    "location": "eigenproblems/lobpcg/#Implementation-Details-1",
    "page": "LOBPCG",
    "title": "Implementation Details",
    "category": "section",
    "text": "A LOBPCGIterator is created to pre-allocate all the memory required by the method using the constructor LOBPCGIterator(A, B, largest, X, P, C) where A and B are the matrices from the generalized eigenvalue problem, largest indicates if the problem is a maximum or minimum eigenvalue problem, X is the initial eigenbasis, randomly sampled if not input, where size(X, 2) is the block size bs. P is the preconditioner, nothing by default, and C is the constraints matrix. The desired k eigenvalues are found bs at a time."
},

{
    "location": "eigenproblems/lobpcg/#References-1",
    "page": "LOBPCG",
    "title": "References",
    "category": "section",
    "text": "Implementation is based on [Knyazev1993] and [Scipy].[Knyazev1993]: Andrew V. Knyazev. \"Toward the Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method\" SIAM Journal on Scientific Computing, 23(2):517–541 2001.[Scipy]: See Scipy LOBPCG implementation"
},

{
    "location": "svd/svdl/#",
    "page": "SVDL",
    "title": "SVDL",
    "category": "page",
    "text": ""
},

{
    "location": "svd/svdl/#SVDL-1",
    "page": "SVDL",
    "title": "Golub-Kahan-Lanczos (SVDL)",
    "category": "section",
    "text": "The SVDL method computes a partial, approximate SVD decomposition of a general linear operator A."
},

{
    "location": "svd/svdl/#IterativeSolvers.svdl",
    "page": "SVDL",
    "title": "IterativeSolvers.svdl",
    "category": "function",
    "text": "svdl(A) -> Σ, L, [history]\n\nCompute some singular values (and optionally vectors) using Golub-Kahan-Lanczos bidiagonalization [Golub1965] with thick restarting [Wu2000].\n\nIf log is set to true is given, method will output a tuple X, L, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return X, L.\n\nArguments\n\nA : The matrix or matrix-like object whose singular values are desired.\n\nKeywords\n\nnsv::Int = 6: number of singular values requested;\nv0 = random unit vector: starting guess vector in the domain of A. The length of q should be the number of columns in A;\nk::Int = 2nsv: maximum number of Lanczos vectors to compute before restarting;\nj::Int = nsv: number of vectors to keep at the end of the restart. We don\'t recommend j < nsv;\nmaxiter::Int = minimum(size(A)): maximum number of iterations to run;\nverbose::Bool = false: print information at each iteration;\ntol::Real = √eps(): maximum absolute error in each desired singular value;\nreltol::Real=√eps(): maximum error in each desired singular value relative to the estimated norm of the input matrix;\nmethod::Symbol=:ritz: restarting algorithm to use. Valid choices are:\n:ritz: Thick restart with Ritz values [Wu2000].\n:harmonic: Restart with harmonic Ritz values [Baglama2005].\nvecs::Symbol = :none: singular vectors to return.\n:both: Both left and right singular vectors are returned.\n:left: Only the left singular vectors are returned.\n:right: Only the right singular vectors are returned.\n:none: No singular vectors are returned.\ndolock::Bool=false: If true, locks converged Ritz values, removing them from the Krylov subspace being searched in the next macroiteration;\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nReturn values\n\nif log is false\n\nΣ: list of the desired singular values if vecs == :none (the default), otherwise returns an SVD object with the desired singular vectors filled in;\nL: computed partial factorizations of A.\n\nif log is true\n\nΣ: list of the desired singular values if vecs == :none (the default),\n\notherwise returns an SVD object with the desired singular vectors filled in;\n\nL: computed partial factorizations of A;\nhistory: convergence history.\n\nConvergenceHistory keys\n\n:betas => betas: The history of the computed betas.\n:Bs => Bs: The history of the computed projected matrices.\n:ritz => ritzvalhist: Ritz values computed at each iteration.\n:conv => convhist: Convergence data.\n\n[Golub1965]: Golub, Gene, and William Kahan. \"Calculating the singular values and pseudo-inverse of a matrix.\" Journal of the Society for Industrial and Applied Mathematics, Series B: Numerical Analysis 2.2 (1965): 205-224.\n\n[Wu2000]: Wu, Kesheng, and Horst Simon. \"Thick-restart Lanczos method for large symmetric eigenvalue problems.\" SIAM Journal on Matrix Analysis and Applications 22.2 (2000): 602-616.\n\n\n\n\n\n"
},

{
    "location": "svd/svdl/#Usage-1",
    "page": "SVDL",
    "title": "Usage",
    "category": "section",
    "text": "svdl"
},

{
    "location": "svd/svdl/#Implementation-details-1",
    "page": "SVDL",
    "title": "Implementation details",
    "category": "section",
    "text": "The implementation of thick restarting follows closely that of SLEPc as described in [Hernandez2008]. Thick restarting can be turned off by setting k = maxiter, but most of the time this is not desirable.The singular vectors are computed directly by forming the Ritz vectors from the product of the Lanczos vectors L.P/L.Q and the singular vectors of L.B. Additional accuracy in the singular triples can be obtained using inverse iteration.[Hernandez2008]: Vicente Hernández, José E. Román, and Andrés Tomás. \"A robust and efficient parallel SVD solver based on restarted Lanczos bidiagonalization.\" Electronic Transactions on Numerical Analysis 31 (2008): 68-85."
},

{
    "location": "iterators/#",
    "page": "The iterator approach",
    "title": "The iterator approach",
    "category": "page",
    "text": ""
},

{
    "location": "iterators/#Iterators-1",
    "page": "The iterator approach",
    "title": "Iterative solvers as iterators",
    "category": "section",
    "text": "In advanced use cases you might want to access the internal data structures of the solver, inject code to be run after each iteration, have total control over allocations or reduce overhead in initialization. The iterator approach of IterativeSolvers.jl makes this possible.note: Note\nAt this point BiCGStab(l), CG, Chebyshev, GMRES, MINRES and the stationary methods are implemented as iterators. However, the package does not yet export the iterators and helper methods themselves."
},

{
    "location": "iterators/#How-iterators-are-implemented-1",
    "page": "The iterator approach",
    "title": "How iterators are implemented",
    "category": "section",
    "text": "The solvers listed above are basically a thin wrapper around an iterator. Among other things, they initialize the iterable, loop through the iterator and return the result:function my_solver!(x, A, b)\n    iterable = MySolverIterable(x, A, b)\n    for item in iterable end\n    return iterable.x\nendRather than calling my_solver!(x, A, b), you could also initialize the iterable yourself and perform the for loop."
},

{
    "location": "iterators/#Example:-avoiding-unnecessary-initialization-1",
    "page": "The iterator approach",
    "title": "Example: avoiding unnecessary initialization",
    "category": "section",
    "text": "The Jacobi method for SparseMatrixCSC has some overhead in intialization; not only do we need to allocate a temporary vector, we also have to search for indices of the diagonal (and check their values are nonzero). The current implementation initializes the iterable as:jacobi_iterable(x, A::SparseMatrixCSC, b; maxiter::Int = 10) =\n    JacobiIterable{eltype(x), typeof(x)}(OffDiagonal(A, DiagonalIndices(A)), x, similar(x), b, maxiter)Now if you apply Jacobi iteration multiple times with the same matrix for just a few iterations, it makes sense to initialize the iterable only once and reuse it afterwards:A = sprand(10_000, 10_000, 10 / 10_000) + 20I\nb1 = rand(10_000)\nb2 = rand(10_000)\nx = rand(10_000)\n\nmy_iterable = IterativeSolvers.jacobi_iterable(x, A, b1, maxiter = 4)\n\nfor item in my_iterable \n    println(\"Iteration for rhs 1\")\nend\n\n@show norm(b1 - A * x) / norm(b1)\n\n# Copy the next right-hand side into the iterable\ncopyto!(my_iterable.b, b2)\n\nfor item in my_iterable\n    println(\"Iteration for rhs 2\")\nend\n\n@show norm(b2 - A * x) / norm(b2)This would output:Iteration for rhs 1\nIteration for rhs 1\nIteration for rhs 1\nIteration for rhs 1\nnorm(b1 - A * x) / norm(b1) = 0.08388528015119746\nIteration for rhs 2\nIteration for rhs 2\nIteration for rhs 2\nIteration for rhs 2\nnorm(b2 - A * x) / norm(b2) = 0.0003681972775644809"
},

{
    "location": "iterators/#Other-use-cases-1",
    "page": "The iterator approach",
    "title": "Other use cases",
    "category": "section",
    "text": "Other use cases include: computing the (harmonic) Ritz values from the Hessenberg matrix in GMRES;\ncomparing the approximate residual of methods such as GMRES and BiCGStab(l) with the true residual during the iterations;\nupdating a preconditioner in flexible methods."
},

{
    "location": "about/CONTRIBUTING/#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "about/CONTRIBUTING/#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "Contributions are always welcome, as are feature requests and suggestions. Feel free to open an issue or pull request at any time.It is important to note that almost every method in the package has documentation, to know what it does simply use ?<method> in the terminal.julia> using IterativeSolvers\n\nhelp?> IterativeSolvers.Adivtype\n  Adivtype(A, b)\n\n  Determine type of the division of an element of b against an element of A:\n\n  typeof(one(eltype(b))/one(eltype(A)))"
},

{
    "location": "about/CONTRIBUTING/#Setting-workspace-up-1",
    "page": "Contributing",
    "title": "Setting workspace up",
    "category": "section",
    "text": "Julia\'s internal package manager makes it easy to install and modify packages from Github. Any package hosted on Github can be installed via Pkg.clone by providing the repository\'s URL, so installing a fork on your system is a simple task.Pkg.clone(\"https://github.com/johndoe/IterativeSolvers.jl\")It is to note here if you have the original package installed the fork will replace it, this is not a problem.Now find your fork\'s location.Pkg.dir(\"IterativeSolvers\")Once there you will notice you are on the master branch, whenever a package is imported Julia will use the code in the current branch, this means checking out other git branches will let you use/test whatever there is."
},

{
    "location": "about/CONTRIBUTING/#Adding-or-modifying-iterative-methods-1",
    "page": "Contributing",
    "title": "Adding or modifying iterative methods",
    "category": "section",
    "text": "Each iterative method method must log information using the inner ConvergenceHistory type.There are two types of ConvergenceHistory: plain and restarted. The only real difference between the two is how they are plotted and how the number of restarts is calculated, everything else is the same.Before logging information space must always be reserved.log = ConvergenceHistory()\nlog[:tol] = tol\nreserve!(log,:betas, maxiter) # Vector of length maxiter\nreserve!(log,:conv, maxiter, T=BitArray) # Vector of length maxiter\nreserve!(log,:ritz, maxiter, k) # Matrix of size (maxiter, k)To store information at each iteration use push!.push!(log, :conv, conv)\npush!(log, :ritz, F[:S][1:k])\npush!(log, :betas, L.β)To advance the log index to the next iteration use nextiter!.nextiter!(log)A more detailed explanation of all the functions is in both the public and internal documentation of ConvergenceHistory.The most rich example of the usage of ConvergenceHistory is in svdl."
},

{
    "location": "about/license/#",
    "page": "License",
    "title": "License",
    "category": "page",
    "text": ""
},

{
    "location": "about/license/#License-(MIT)-1",
    "page": "License",
    "title": "License (MIT)",
    "category": "section",
    "text": "Copyright (c) 2013--2016 The Julia Language\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of\nthis software and associated documentation files (the \"Software\"), to deal in\nthe Software without restriction, including without limitation the rights to\nuse, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of\nthe Software, and to permit persons to whom the Software is furnished to do so,\nsubject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS\nFOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR\nCOPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER\nIN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN\nCONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

]}
