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
    "text": "IterativeSolvers.jl is a Julia package that provides iterative algorithms for solving linear systems, eigensystems, and singular value problems. The purpose of this package is to provide efficient Julia implementations for iterative methods. The package aims to accept a wide variety of input types and that's why most arguments don't specify a specific type, however this is still in progress.For bug reports, feature requests and questions please submit an issue. If you're interested in contributing, please see the Contributing guide.For more information on future methods have a look at the package roadmap on deterministic methods, for randomized algorithms check here."
},

{
    "location": "index.html#Linear-Solvers-1",
    "page": "Home",
    "title": "Linear Solvers",
    "category": "section",
    "text": "Stationary methodsJacobi\nGauss-Seidel\nSuccessive over relaxation\nSymmetric successive over relaxationNon stationary methodsIDRS\nLSMR\nLSQR\nConjugate gradients\nChebyshev iteration\nGeneralized minimal residual method (with restarts)"
},

{
    "location": "index.html#Eigen-Solvers-1",
    "page": "Home",
    "title": "Eigen Solvers",
    "category": "section",
    "text": "Simple eigenpair iterationsPower iteration\nInverse power iterationHermitianLanczosSimple Lanczos"
},

{
    "location": "index.html#Singular-Value-Decomposition-1",
    "page": "Home",
    "title": "Singular Value Decomposition",
    "category": "section",
    "text": "Golub-Kahan-Lanczos"
},

{
    "location": "index.html#Randomized-1",
    "page": "Home",
    "title": "Randomized",
    "category": "section",
    "text": "Condition number estimate\nExtremal eigenvalue estimates\nNorm estimate\nRandomized singular value decomposition"
},

{
    "location": "index.html#Documentation-Outline-1",
    "page": "Home",
    "title": "Documentation Outline",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\n    \"user_manual.md\",\n]\nDepth = 2"
},

{
    "location": "index.html#Library-1",
    "page": "Home",
    "title": "Library",
    "category": "section",
    "text": "Pages = [\"library/public.md\", \"library/internal.md\"]\nDepth = 2"
},

{
    "location": "index.html#main-index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Functions-1",
    "page": "Home",
    "title": "Functions",
    "category": "section",
    "text": "Pages = [\"library/public.md\", \"library/internals.md\"]\nOrder = [:function]"
},

{
    "location": "index.html#Types-1",
    "page": "Home",
    "title": "Types",
    "category": "section",
    "text": "Pages = [\"library/public.md\", \"library/internals.md\"]\nOrder = [:type]"
},

{
    "location": "user_manual.html#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "user_manual.html#Manual-1",
    "page": "Manual",
    "title": "Manual",
    "category": "section",
    "text": ""
},

{
    "location": "user_manual.html#Installation-1",
    "page": "Manual",
    "title": "Installation",
    "category": "section",
    "text": "The package can be installed with a simple instruction.julia> Pkg.add(\"IterativeSolvers\")After installing the package, if you wish to use the latest features of the package you must switch to the master branch with Pkg.clone.julia> Pkg.checkout(\"IterativeSolvers\")"
},

{
    "location": "user_manual.html#Interface-1",
    "page": "Manual",
    "title": "Interface",
    "category": "section",
    "text": "All linear-algebraic routines will take as input a linear operator A that maps vectors to vectors. A is not explicitly typed, but must either be a KrylovSubspace or support multiplication * or function composition (apply) that behave as necessary to produce the correct mapping on the vector space.A custom type for A may be specified. The following interface is expected to be defined on A:A*v is defined and computes the matrix-vector product on a v::Vector.eltype(A) is defined and returns the element type implicit in the equivalent matrix representation of A.size(A, d) is defined and returns the nominal dimensions along the dth axis in the equivalent matrix representation of A."
},

{
    "location": "user_manual.html#Solvers-1",
    "page": "Manual",
    "title": "Solvers",
    "category": "section",
    "text": "All linear solvers have a common function declaration (with a few exceptions).solver(A, b::Vector; kwargs...)\nsolver!(x, A, b::Vector; kwargs...)In the case of eigenproblems or singular value decompositions:eigsolver(A; kwargs...)\neigsolver!(x, A; kwargs...)A is a linear operator as described above.b is the vector to be solved.x is a vector for the initial guess. In the case of a mutating call this parameter will be overwritten. (Mutating functions end with !)Output will be the solution to the system."
},

{
    "location": "user_manual.html#Additional-arguments-1",
    "page": "Manual",
    "title": "Additional arguments",
    "category": "section",
    "text": "Keyword names will vary depending on the method, however some of them will always have the same spelling:tol: stopping tolerance of the method. When a method accepts more than one tolerance they are enumerated  with a letter prefix, e.g atol, btol, ctol, etc.verbose: print information about the running method.maxiter: maximum number of allowed iterations.Pl: left preconditioner. (When applicable)Pr: right preconditioner. (When applicable)log::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.plot: plot information relevant to the method. (Only for Master version)"
},

{
    "location": "user_manual.html#log-keyword-1",
    "page": "Manual",
    "title": "log keyword",
    "category": "section",
    "text": "All solvers contain the log keyword. This is to be used when obtaining more information is required, to use it place the set log to true.x, ch = cg(Master, rand(10, 10), rand(10) log=true)\nsvd, L, ch = svdl(Master, rand(100, 100), plot=true, log=true)The function will now return one more parameter of type ConvergenceHistory.Note:  Keyword argument plot is only available when log is set."
},

{
    "location": "user_manual.html#ConvergenceHistory-1",
    "page": "Manual",
    "title": "ConvergenceHistory",
    "category": "section",
    "text": "A ConvergenceHistory instance stores information of a solver.Number of iterations.ch.itersConvergence status.ch.isconvergedStopping tolerances. (A Symbol key is needed to access)ch[:tol]Maximum number of iterations per restart. (Only on restarted methods)nrests(ch)Number of matrix-vectors and matrix-transposed-vector products.nprods(ch)Data stored on each iteration, accessed information can be either a vector or matrix. This data can be a lot of things, most commonly residual. (A Symbol key is needed to access)ch[:resnorm] #Vector or Matrix\nch[:resnorm, x] #Vector or Matrix element\nch[:resnorm, x, y] #Matrix elementThe available keys of each method is described in the Public Documentation."
},

{
    "location": "user_manual.html#Plotting-1",
    "page": "Manual",
    "title": "Plotting",
    "category": "section",
    "text": "ConvergeHistory provides a recipe to use with the package Plots.jl, this makes it really easy to plot on different plot backends. There are two recipes provided:One for the whole ConvergenceHistory.plot(ch)The other one to plot data binded to a key._, ch = gmres(rand(10,10), rand(10), maxiter = 100, log=true)\nplot(ch, :resnorm, sep = :blue)Plot additional keywordssep::Symbol = :white: color of the line separator in restarted methods."
},

{
    "location": "user_manual.html#KrylovSubspace-1",
    "page": "Manual",
    "title": "KrylovSubspace",
    "category": "section",
    "text": "When KrylovSubspace is supported by the method, A can be replaced by an instance K of it. To check if a certain function supports it use   methods or ? to find out, if there is a K in the argument list then it does.julia> ?cg!\n\n    cg!(x, A, b)\n    cg!(x, K, b)\n\n    ...KrylovSubspace is only allowed on functions that accept x as an argument, that's why cg doesn't allow it. This is also true when x is a keyword as is the case of powm.julia> ?cg!\n\n    powm(A)\n    powm(K)\n\n    ...The KrylovSubspace type collects information on the Krylov subspace generated over the course of an iterative Krylov solver.Recall that the Krylov subspace of order r given a starting vector b and a linear operator A is spanned by the vectors [b, A*b, A^2*b,... A^(r-1)*b]. Many modern iterative solvers work on Krylov spaces which expand until they span enough of the range of A for the solution vector to converge. Most practical algorithms, however, will truncate the order of the Krylov subspace to save memory and control the accumulation of roundoff errors. Additionally, they do not work directly with the raw iterates A^n*b, but will orthogonalize subsequent vectors as they are computed so as to improve numerical stability. KrylovSubspaces provide a natural framework to express operations such as a (possibly non-orthonormalized) basis for the Krylov subspace, retrieving the next vector in the subspace, and orthogonalizing an arbitrary vector against (part or all of) the subspace.The implementation of KrylovSubspace in this package differs from standard textbook presentations of iterative solvers. First, the KrylovSubspace type shows clearly the relationship between the linear operator A and the sequence of basis vectors for the Krylov subspace that is generated by each new iteration. Second, the grouping together of basis vectors also makes the orthogonalization steps in each iteration routine clearer. Third, unlike in other languages, the powerful type system of Julia allows for a simplified notation for iterative algorithms without compromising on performance, and enhances code reuse across multiple solvers."
},

{
    "location": "user_manual.html#Constructors-1",
    "page": "Manual",
    "title": "Constructors",
    "category": "section",
    "text": "A KrylovSubspace can be initialized by three constructors, depending on the type of the linear operator:KrylovSubspace{T}(A::AbstractMatrix{T}, [order::Int, v::Vector{Vector{T}}])\n\nKrylovSubspace{T}(A::KrylovSubspace{T}, [order::Int, v::Vector{Vector{T}}])\n\nKrylovSubspace(A, n::Int, [order::Int, v::Vector{Vector{T}}])A: The linear operator associated with the KrylovSubspace.order: the order of the KrylovSubspace, i.e. the maximal number of Krylov vectors to remember.n: the dimensionality of the underlying vector space T^n.v: the iterable collection of Krylov vectors (of maximal length order).The dimensionality of the underlying vector space is automatically inferred where possible, i.e. when the linear operator is an AbstractMatrix or KrylovSubspace.Note: second constructor destroys the old KrylovSubspace."
},

{
    "location": "user_manual.html#Orthogonalization-1",
    "page": "Manual",
    "title": "Orthogonalization",
    "category": "section",
    "text": "Orthogonalizing the basis vectors for a KrylovSubspace is crucial for numerical stability, and is a core low-level operation for many iterative solvers.orthogonalize{T}(v::Vector{T}, K::KrylovSubspace{T}, [p::Int]; [method::Symbol], [normalize::Bool])v: vector to orthogonalize.K: KrylovSubspace to orthogonalize against.p: number of Krylov vectors to orthogonalize against. (Default is all available)method: Orthogonalization method. Currently supported methods are::GramSchmidt\n:ModifiedGramSchmidt (default)\n:Householder."
},

{
    "location": "user_manual.html#Define-function-as-matrix-1",
    "page": "Manual",
    "title": "Define function as matrix",
    "category": "section",
    "text": "Suppose you have a function implementing multiplication of b by A. This function can have any syntax, but for the purposes of illustration let's suppose it's defined as:mulbyA!(output, b, Adata)Where Adata might be some parameters that your function needs."
},

{
    "location": "library/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "library/public.html#Public-Documentation-1",
    "page": "Public",
    "title": "Public Documentation",
    "category": "section",
    "text": "Documentation for IterativeSolvers.jl's public interface.Pages = [\"public.md\"]\nDepth = 4"
},

{
    "location": "library/public.html#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "library/public.html#Types-1",
    "page": "Public",
    "title": "Types",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.ConvergenceHistory",
    "page": "Public",
    "title": "IterativeSolvers.ConvergenceHistory",
    "category": "Type",
    "text": "Store general and in-depth information about an iterative method.\n\nFields\n\nmvps::Int: number of matrix vector products.\n\nmtvps::Int: number of transposed matrix-vector products\n\niters::Int: iterations taken by the method.\n\nrestart::T: restart relevant information.\n\nT == Int: iterations per restart.\nT == Void: methods without restarts.\n\nisconverged::Bool: convergence of the method.\n\ndata::Dict{Symbol,Any}: Stores all the information stored during the method execution. It stores tolerances, residuals and other information, e.g. ritz values in svdl.\n\nConstructors\n\nConvergenceHistory()\nConvergenceHistory(restart)\n\nCreate ConvergenceHistory with empty fields.\n\nArguments\n\nrestart: number of iterations per restart.\n\nPlots\n\nSupports plots using the Plots.jl package via a type recipe. Vectors are ploted as series and matrices as scatterplots.\n\nImplements\n\nBase: getindex, setindex!, push!\n\n\n\n"
},

{
    "location": "library/public.html#ConvergenceHistory-1",
    "page": "Public",
    "title": "ConvergenceHistory",
    "category": "section",
    "text": "ConvergenceHistory"
},

{
    "location": "library/public.html#IterativeSolvers.KrylovSubspace",
    "page": "Public",
    "title": "IterativeSolvers.KrylovSubspace",
    "category": "Type",
    "text": "Collection of information on the Krylov subspace generated over the course of an iterative Krylov solver.\n\nFields\n\nA: linear operator that generating the subspace.\n\nn::Int: dimension of problem.\n\norder::Int: order of subspace (maximum size).\n\nv::Vector{Vector}: Krylov vectors.\n\nmvps::Int: count of matrix-vector products.\n\nConstructors\n\nKrylovSubspace(A, order)\nKrylovSubspace(A, n, order)\nKrylovSubspace(A, n, order, T)\nKrylovSubspace(A, n, order, v)\n\nCreate new KrylovSubspace.\n\nArguments\n\nT::Type: type of the elements inside the Krylov vectors (eltype(v)).\n\nv::Vector{Vector}: initial Krylov subspace.\n\nKrylovSubspace(K, order=size(A,2), v=[])\n\nReset given KrylovSubspace A.\n\nArguments\n\nK::KrylovSubspace: Krylov subspace to reset.\n\nImplements\n\nBase: append!, eye, size.\n\n\n\n"
},

{
    "location": "library/public.html#KrylovSubspace-1",
    "page": "Public",
    "title": "KrylovSubspace",
    "category": "section",
    "text": "KrylovSubspace"
},

{
    "location": "library/public.html#Linear-Solvers-1",
    "page": "Public",
    "title": "Linear Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.jacobi",
    "page": "Public",
    "title": "IterativeSolvers.jacobi",
    "category": "Function",
    "text": "jacobi(A, b)\n\nSolve A*x=b with the Jacobi method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.jacobi!",
    "page": "Public",
    "title": "IterativeSolvers.jacobi!",
    "category": "Function",
    "text": "jacobi!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b with the Jacobi method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Jacobi-1",
    "page": "Public",
    "title": "Jacobi",
    "category": "section",
    "text": "jacobi\njacobi!"
},

{
    "location": "library/public.html#IterativeSolvers.gauss_seidel",
    "page": "Public",
    "title": "IterativeSolvers.gauss_seidel",
    "category": "Function",
    "text": "gauss_seidel(A, b)\n\nSolve A*x=b with the Gauss-Seidel method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.gauss_seidel!",
    "page": "Public",
    "title": "IterativeSolvers.gauss_seidel!",
    "category": "Function",
    "text": "gauss_seidel!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b with the Gauss-Seidel method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Gauss-Seidel-1",
    "page": "Public",
    "title": "Gauss-Seidel",
    "category": "section",
    "text": "gauss_seidel\ngauss_seidel!"
},

{
    "location": "library/public.html#IterativeSolvers.sor",
    "page": "Public",
    "title": "IterativeSolvers.sor",
    "category": "Function",
    "text": "sor(A, b, ω)\n\nSolve A*x=b with the successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.sor!",
    "page": "Public",
    "title": "IterativeSolvers.sor!",
    "category": "Function",
    "text": "sor!(x, A, b, ω)\n\nOverwrite x.\n\nSolve A*x=b with the successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Successive-over-relaxation-1",
    "page": "Public",
    "title": "Successive over relaxation",
    "category": "section",
    "text": "sor\nsor!"
},

{
    "location": "library/public.html#IterativeSolvers.ssor",
    "page": "Public",
    "title": "IterativeSolvers.ssor",
    "category": "Function",
    "text": "ssor(A, b, ω)\n\nSolve A*x=b with the symmetric successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.ssor!",
    "page": "Public",
    "title": "IterativeSolvers.ssor!",
    "category": "Function",
    "text": "ssor!(x, A, b, ω)\n\nOverwrite x.\n\nSolve A*x=b with the symmetric successive overrelaxation method.\n\nch is a ConvergenceHistory object. Otherwise it will only return x. If log is set to true is given, method will output a tuple x, ch. Where The plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nshift::Number=0: shift to be applied to matrix A.\n\nA: linear operator.\n\nKeywords\n\ntol::Real = size(A,2)^3*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^2: maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Symmetric-successive-over-relaxation-1",
    "page": "Public",
    "title": "Symmetric successive over relaxation",
    "category": "section",
    "text": "ssor\nssor!"
},

{
    "location": "library/public.html#IterativeSolvers.idrs",
    "page": "Public",
    "title": "IterativeSolvers.idrs",
    "category": "Function",
    "text": "idrs(A, b)\n\nSolve A*x=b using the induced dimension reduction method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\nPr = 1: left preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nrestart::Integer = min(20,length(b)): maximum number of iterations per restart.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.idrs!",
    "page": "Public",
    "title": "IterativeSolvers.idrs!",
    "category": "Function",
    "text": "idrs!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b using the induced dimension reduction method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\nPr = 1: left preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nrestart::Integer = min(20,length(b)): maximum number of iterations per restart.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IDRS-1",
    "page": "Public",
    "title": "IDRS",
    "category": "section",
    "text": "idrs\nidrs!References[1] IDR(s): a family of simple and fast algorithms for solving large\n    nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen\n    SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, 2008\n[2] Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits\n    Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld\n    ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011\n[3] This file is a translation of the following MATLAB implementation:\n    http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m\n[4] IDR(s)' webpage http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html"
},

{
    "location": "library/public.html#IterativeSolvers.lsmr",
    "page": "Public",
    "title": "IterativeSolvers.lsmr",
    "category": "Function",
    "text": "lsmr(A, b)\n\nMinimize ||Ax-b||^2 + λ^2 ||x||^2 for A*x=b.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying MINRES to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.lsmr!",
    "page": "Public",
    "title": "IterativeSolvers.lsmr!",
    "category": "Function",
    "text": "lsmr!(x, A, b)\n\nOverwrite x.\n\nMinimize ||Ax-b||^2 + λ^2 ||x||^2 for A*x=b.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying MINRES to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#LSMR-1",
    "page": "Public",
    "title": "LSMR",
    "category": "section",
    "text": "lsmr\nlsmr!ReferencesAdapted from: http://web.stanford.edu/group/SOL/software/lsmr/"
},

{
    "location": "library/public.html#IterativeSolvers.lsqr",
    "page": "Public",
    "title": "IterativeSolvers.lsqr",
    "category": "Function",
    "text": "lsqr(A, b)\n\nLSQR solves Ax = b or min ||b - Ax||^2 if damp = 0, or   min ||(b) - (  A   )x||   otherwise.          ||(0)   (damp*I) ||^2.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying CG to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.lsqr!",
    "page": "Public",
    "title": "IterativeSolvers.lsqr!",
    "category": "Function",
    "text": "lsqr!(x, A, b)\n\nOverwrite x.\n\nLSQR solves Ax = b or min ||b - Ax||^2 if damp = 0, or   min ||(b) - (  A   )x||   otherwise.          ||(0)   (damp*I) ||^2.\n\nThe method is based on the Golub-Kahan bidiagonalization process. It is algebraically equivalent to applying CG to the normal equation (ATA+λ2I)x=ATb, but has better numerical properties, especially if A is ill-conditioned.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nλ::Number = 0: lambda.\n\natol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).\n\nconlim::Number = 1e8: stopping tolerance.  lsmr terminates if an estimate of cond(A) exceeds conlim.  For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say).  For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:atol => ::Real: atol stopping tolerance.\n\n:btol => ::Real: btol stopping tolerance.\n\n:ctol => ::Real: ctol stopping tolerance.\n\n:anorm => ::Real: anorm.\n\n:rnorm => ::Real: rnorm.\n\n:cnorm => ::Real: cnorm.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#LSQR-1",
    "page": "Public",
    "title": "LSQR",
    "category": "section",
    "text": "lsqr\nlsqr!ReferencesAdapted from: http://web.stanford.edu/group/SOL/software/lsqr/\n\n1. C. C. Paige and M. A. Saunders (1982a).\n    LSQR: An algorithm for sparse linear equations and sparse least squares,\n    ACM TOMS 8(1), 43-71.\n\n2. C. C. Paige and M. A. Saunders (1982b).\n    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,\n    ACM TOMS 8(2), 195-209.\n\n3. M. A. Saunders (1995).  Solution of sparse rectangular systems using\n    LSQR and CRAIG, BIT 35, 588-604."
},

{
    "location": "library/public.html#IterativeSolvers.cg",
    "page": "Public",
    "title": "IterativeSolvers.cg",
    "category": "Function",
    "text": "cg(A, b)\n\nSolve A*x=b with the conjugate gradients method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x. The plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\ntol::Real = size(A,2)*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.cg!",
    "page": "Public",
    "title": "IterativeSolvers.cg!",
    "category": "Function",
    "text": "cg!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b with the conjugate gradients method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x. The plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\ntol::Real = size(A,2)*eps(): stopping tolerance.\n\nmaxiter::Integer = size(A,2): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Conjugate-gradients-1",
    "page": "Public",
    "title": "Conjugate gradients",
    "category": "section",
    "text": "cg\ncg!"
},

{
    "location": "library/public.html#IterativeSolvers.chebyshev",
    "page": "Public",
    "title": "IterativeSolvers.chebyshev",
    "category": "Function",
    "text": "chebyshev!(A, b, λmin, λmax)\n\nSolve A*x=b using the chebyshev method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPr = 1: right preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^3: maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only with Master version)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.chebyshev!",
    "page": "Public",
    "title": "IterativeSolvers.chebyshev!",
    "category": "Function",
    "text": "chebyshev!(x, A, b, λmin, λmax)\n\nOverwrite x.\n\nSolve A*x=b using the chebyshev method.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPr = 1: right preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nmaxiter::Integer = size(A,2)^3: maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only with Master version)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Chebyshev-iteration-1",
    "page": "Public",
    "title": "Chebyshev iteration",
    "category": "section",
    "text": "chebyshev\nchebyshev!"
},

{
    "location": "library/public.html#IterativeSolvers.gmres",
    "page": "Public",
    "title": "IterativeSolvers.gmres",
    "category": "Function",
    "text": "gmres(A, b)\n\nSolve A*x=b using the generalized minimal residual method with restarts.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\nPr = 1: left preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nrestart::Integer = min(20,length(b)): maximum number of iterations per restart.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance. :resnom => ::Vector: residual norm at each iteration.\n\nReferences\n\nhttp://www.netlib.org/templates/templates.pdf 2.3.4 Generalized Minimal Residual (GMRES)\n\nhttp://www.netlib.org/lapack/lawnspdf/lawn148.pdf Givens rotation based on Algorithm 1\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.gmres!",
    "page": "Public",
    "title": "IterativeSolvers.gmres!",
    "category": "Function",
    "text": "gmres!(x, A, b)\n\nOverwrite x.\n\nSolve A*x=b using the generalized minimal residual method with restarts.\n\nIf log is set to true is given, method will output a tuple x, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return x.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nx: initial guess, overwrite final estimation.\n\nA: linear operator.\n\nb: right hand side.\n\nKeywords\n\nPl = 1: left preconditioner of the method.\n\nPr = 1: left preconditioner of the method.\n\ntol::Real = sqrt(eps()): stopping tolerance.\n\nrestart::Integer = min(20,length(b)): maximum number of iterations per restart.\n\nmaxiter::Integer = min(20,length(b)): maximum number of iterations.\n\nverbose::Bool = false: print method information.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nx: approximated solution.\n\nif log is true\n\nx: approximated solution.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance. :resnom => ::Vector: residual norm at each iteration.\n\nReferences\n\nhttp://www.netlib.org/templates/templates.pdf 2.3.4 Generalized Minimal Residual (GMRES)\n\nhttp://www.netlib.org/lapack/lawnspdf/lawn148.pdf Givens rotation based on Algorithm 1\n\n\n\n"
},

{
    "location": "library/public.html#Generalized-minimal-residual-method-(with-restarts)-1",
    "page": "Public",
    "title": "Generalized minimal residual method (with restarts)",
    "category": "section",
    "text": "gmres\ngmres!"
},

{
    "location": "library/public.html#Eigen-Solvers-1",
    "page": "Public",
    "title": "Eigen Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.powm",
    "page": "Public",
    "title": "IterativeSolvers.powm",
    "category": "Function",
    "text": "powm(A)\n\nFind biggest eigenvalue of A and its associated eigenvector using the power method.\n\nIf log is set to true is given, method will output a tuple eig, v, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return eig, v.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nK::KrylovSubspace: krylov subspace.\n\nA: linear operator.\n\nKeywords\n\nx = random unit vector: initial eigenvector guess.\n\ntol::Real = eps()*size(A,2)^3: stopping tolerance.\n\nmaxiter::Integer = size(A,2): maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\neig::Real: eigen value\n\nv::Vector: eigen vector\n\nif log is true\n\neig::Real: eigen value\n\nv::Vector: eigen vector\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Power-iteration-1",
    "page": "Public",
    "title": "Power iteration",
    "category": "section",
    "text": "powm"
},

{
    "location": "library/public.html#IterativeSolvers.invpowm",
    "page": "Public",
    "title": "IterativeSolvers.invpowm",
    "category": "Function",
    "text": "invpowm(A)\n\nFind closest eigenvalue of A to shift and its associated eigenvector using the inverse power iteration method.\n\nIf log is set to true is given, method will output a tuple eig, v, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return eig, v.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nK::KrylovSubspace: krylov subspace.\n\nA: linear operator.\n\nKeywords\n\nshift::Number=0: shift to be applied to matrix A.\n\nx = random unit vector: initial eigenvector guess.\n\ntol::Real = eps()*size(A,2)^3: stopping tolerance.\n\nmaxiter::Integer = size(A,2): maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\neig::Real: eigen value\n\nv::Vector: eigen vector\n\nif log is true\n\neig::Real: eigen value\n\nv::Vector: eigen vector\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Inverse-power-iteration-1",
    "page": "Public",
    "title": "Inverse power iteration",
    "category": "section",
    "text": "invpowm"
},

{
    "location": "library/public.html#IterativeSolvers.eiglancz",
    "page": "Public",
    "title": "IterativeSolvers.eiglancz",
    "category": "Function",
    "text": "eiglancz(A)\n\nFind the most useful eigenvalues using the lanczos method.\n\nIf log is set to true is given, method will output a tuple eigs, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return eigs.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA: linear operator.\n\nKeywords\n\nneigs::Int = size(A,1): number of eigen values.\n\ntol::Real = size(A,1)^3*eps(): stopping tolerance.\n\nmaxiter::Integer=size(A,1): maximum number of iterations.\n\nverbose::Bool = false: verbose flag.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\neigs::Vector: eigen values.\n\nif log is true\n\neigs::Vector: eigen values.\n\nch: convergence history.\n\nConvergenceHistory keys\n\n:tol => ::Real: stopping tolerance.\n\n:resnom => ::Vector: residual norm at each iteration.\n\n\n\n"
},

{
    "location": "library/public.html#Simple-Lanczos-1",
    "page": "Public",
    "title": "Simple Lanczos",
    "category": "section",
    "text": "eiglancz"
},

{
    "location": "library/public.html#IterativeSolvers.svdl",
    "page": "Public",
    "title": "IterativeSolvers.svdl",
    "category": "Function",
    "text": "svdl(A)\n\nCompute some singular values (and optionally vectors) using Golub-Kahan-Lanczos bidiagonalization cite{Golub1965} with thick restarting cite{Wu2000}.\n\nIf log is set to true is given, method will output a tuple X, L, ch. Where ch is a ConvergenceHistory object. Otherwise it will only return X, L.\n\nThe plot attribute can only be used when log is set version.\n\nArguments\n\nA : The matrix or matrix-like object whose singular values are desired.\n\nKeywords\n\nnsv::Int = 6: number of singular values requested.\n\nv0 = random unit vector: starting guess vector in the domain of A. The length of q should be the number of columns in A.\n\nk::Int = 2nsv: maximum number of Lanczos vectors to compute before restarting.\n\nj::Int = nsv: number of vectors to keep at the end of the restart. We don't recommend j < nsv.\n\nmaxiter::Int = minimum(size(A)): maximum number of iterations to run.\n\nverbose::Bool = false: print information at each iteration.\n\ntol::Real = √eps(): maximum absolute error in each desired singular value.\n\nreltol::Real=√eps(): maximum error in each desired singular value relative to the estimated norm of the input matrix.\n\nmethod::Symbol=:ritz: restarting algorithm to use. Valid choices are:\n\n:ritz: Thick restart with Ritz values [Wu2000].\n:harmonic: Restart with harmonic Ritz values [Baglama2005].\n\nvecs::Symbol = :none: singular vectors to return.\n\n:both: Both left and right singular vectors are returned.\n:left: Only the left singular vectors are returned.\n:right: Only the right singular vectors are returned.\n:none: No singular vectors are returned.\n\ndolock::Bool=false: If true, locks converged Ritz values, removing them from the Krylov subspace being searched in the next macroiteration.\n\nlog::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.\n\nplot::Bool = false: plot data. (Only when log is set)\n\nOutput\n\nif log is false\n\nΣ: list of the desired singular values if vecs == :none (the default), otherwise returns an SVD object with the desired singular vectors filled in.\n\nL: computed partial factorizations of A.\n\nif log is true\n\nΣ: list of the desired singular values if vecs == :none (the default), otherwise returns an SVD object with the desired singular vectors filled in.\n\nL: computed partial factorizations of A.\n\nch::ConvergenceHistory: convergence history.\n\nConvergenceHistory keys\n\n:betas => betas: The history of the computed betas.\n\n:Bs => Bs: The history of the computed projected matrices.\n\n:ritz => ritzvalhist: Ritz values computed at each iteration.\n\n:conv => convhist: Convergence data.\n\n\n\n"
},

{
    "location": "library/public.html#Golub-Kahan-Lanczos-1",
    "page": "Public",
    "title": "Golub-Kahan-Lanczos",
    "category": "section",
    "text": "svdlImplementation notesThe implementation of thick restarting follows closely that of SLEPc as described in [Hernandez2008]. Thick restarting can be turned off by setting k = maxiter, but most of the time this is not desirable.The singular vectors are computed directly by forming the Ritz vectors from the product of the Lanczos vectors L.P/L.Q and the singular vectors of L.B. Additional accuracy in the singular triples can be obtained using inverse iteration.References@article{Golub1965,\n    author = {Golub, G. and Kahan, W.},\n    doi = {10.1137/0702016},\n    journal = {Journal of the Society for Industrial and Applied Mathematics\n        Series B Numerical Analysis},\n    volume = 2,\n    number = 2,\n    pages = {205--224},\n    title = {Calculating the Singular Values and Pseudo-Inverse of a Matrix},\n    year = 1965\n}\n\n@article{Wu2000,\n    author = {Wu, Kesheng and Simon, Horst},\n    journal = {SIAM Journal on Matrix Analysis and Applications},\n    number = 2,\n    pages = {602--616},\n    title = {Thick-Restart {L}anczos Method for Large Symmetric Eigenvalue Problems},\n    volume = 22,\n    year = 2000\n}\n\n@article{Baglama2005,\n    author = {Baglama, James and Reichel, Lothar},\n    doi = {10.1137/04060593X},\n    journal = {SIAM Journal on Scientific Computing},\n    number = 1,\n    pages = {19--42},\n    title = {Augmented Implicitly Restarted {L}anczos Bidiagonalization Methods},\n    volume = 27,\n    year = 2005\n}\n\n@article{Hernandez2008,\n    author = {Hern\\'{a}ndez, Vicente and Rom\\'{a}n, Jos\\'{e} E and Tom\\'{a}s,\n    Andr\\'{e}s},\n    journal = {Electronic Transactions on Numerical Analysis},\n    pages = {68--85},\n    title = {A Robust and Efficient Parallel {SVD} Solver based on Restarted\n        {L}anczos Bidiagonalization},\n    url = {http://etna.mcs.kent.edu/volumes/2001-2010/vol31/abstract.php?vol=31\\&pages=68-85},\n    volume = 31,\n    year = 2008\n}"
},

{
    "location": "library/public.html#Randomized-1",
    "page": "Public",
    "title": "Randomized",
    "category": "section",
    "text": ""
},

{
    "location": "library/public.html#IterativeSolvers.rcond",
    "page": "Public",
    "title": "IterativeSolvers.rcond",
    "category": "Function",
    "text": "rcond(A, iters=1)\n\nEstimate matrix condition number randomly.\n\nArguments\n\nA: matrix whose condition number to estimate. Must be square and support premultiply (A*⋅) and solve (A⋅).\n\niters::Int = 1: number of power iterations to run.\n\nKeywords\n\np::Real = 0.05: probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains κ(A) with probability 1 - p.\n\nImplementation note\n\ncite{Dixon1983} originally describes this as a computation that can be done by computing the necessary number of power iterations given p and the desired accuracy parameter θ=y/x. However, these bounds were only derived under the assumptions of exact arithmetic. Empirically, iters≥4 has been seen to result in incorrect results in that the computed interval does not contain the true condition number. This implemention therefore makes iters an explicitly user-controllable parameter from which to infer the accuracy parameter and hence the interval containing κ(A).\n\nReferences\n\ncite[Theorem 2]{Dixon1983}\n\n@article{Dixon1983,\n    author = {Dixon, John D},\n    doi = {10.1137/0720053},\n    journal = {SIAM Journal on Numerical Analysis},\n    number = {4},\n    pages = {812--814},\n    title = {Estimating Extremal Eigenvalues and Condition Numbers of\n	Matrices},\n    volume = {20},\n    year = {1983}\n}\n\n\n\n"
},

{
    "location": "library/public.html#Condition-number-estimate-1",
    "page": "Public",
    "title": "Condition number estimate",
    "category": "section",
    "text": "rcond"
},

{
    "location": "library/public.html#IterativeSolvers.reigmin",
    "page": "Public",
    "title": "IterativeSolvers.reigmin",
    "category": "Function",
    "text": "reigmin(A, iters=1)\n\nEstimate minimal eigenvalue randomly.\n\nArguments\n\nA: Matrix whose maximal eigenvalue to estimate. Must be square and support premultiply (A*⋅).\n\niters::Int=1: Number of power iterations to run. (Recommended: iters ≤ 3)\n\nKeywords\n\np::Real=0.05: Probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains the maximal eigenvalue of A with probability 1 - p.\n\nReferences\n\ncite[Corollary of Theorem 1]{Dixon1983}.\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.reigmax",
    "page": "Public",
    "title": "IterativeSolvers.reigmax",
    "category": "Function",
    "text": "reigmax(A, iters=1)\n\nEstimate maximal eigenvalue randomly.\n\nArguments\n\nA: Matrix whose maximal eigenvalue to estimate. Must be square and support premultiply (A*⋅).\n\niters::Int=1: Number of power iterations to run. (Recommended: iters ≤ 3)\n\nKeywords\n\np::Real=0.05: Probability that estimate fails to hold as an upper bound.\n\nOutput\n\nInterval (x, y) which contains the maximal eigenvalue of A with probability 1 - p.\n\nReferences\n\ncite[Corollary of Theorem 1]{Dixon1983}.\n\n\n\n"
},

{
    "location": "library/public.html#Extremal-eigenvalue-estimates-1",
    "page": "Public",
    "title": "Extremal eigenvalue estimates",
    "category": "section",
    "text": "reigmin\nreigmax"
},

{
    "location": "library/public.html#IterativeSolvers.rnorm",
    "page": "Public",
    "title": "IterativeSolvers.rnorm",
    "category": "Function",
    "text": "rnorm(A, mvps)\n\nCompute a probabilistic upper bound on the norm of a matrix A. ‖A‖ ≤ α √(2/π) maxᵢ ‖Aωᵢ‖ with probability p=α^(-mvps).\n\nArguments\n\nA: matrix whose norm to estimate.\n\nmvps::Int: number of matrix-vector products to compute.\n\nKeywords\n\np::Real=0.05: probability of upper bound failing.\n\nOutput\n\nEstimate of ‖A‖.\n\nSee also\n\nsee rnorms for a different estimator that uses premultiplying by both A and A'.\n\nReferences\n\ncite[Lemma 4.1]{Halko2011}\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.rnorms",
    "page": "Public",
    "title": "IterativeSolvers.rnorms",
    "category": "Function",
    "text": "rnorms(A, iters=1)\n\nEstimate matrix norm randomly using A'A.\n\nCompute a probabilistic upper bound on the norm of a matrix A.\n\nρ = √(‖(A'A)ʲω‖/‖(A'A)ʲ⁻¹ω‖)\n\nwhich is an estimate of the spectral norm of A produced by iters steps of the power method starting with normalized ω, is a lower bound on the true norm by a factor\n\nρ ≤ α ‖A‖\n\nwith probability greater than 1 - p, where p = 4sqrt(n/(iters-1)) α^(-2iters).\n\nArguments\n\nA: matrix whose norm to estimate.\n\niters::Int = 1: mumber of power iterations to perform.\n\nKeywords\n\np::Real = 0.05: probability of upper bound failing.\n\nAt = A': Transpose of A.\n\nOutput\n\nEstimate of ‖A‖.\n\nSee also\n\nsee rnorm for a different estimator that does not require premultiplying by A'\n\nReferences\n\nAppendix of cite{Liberty2007}.\n\n@article{Liberty2007,\n    authors = {Edo Liberty and Franco Woolfe and Per-Gunnar Martinsson\n    and Vladimir Rokhlin and Mark Tygert},\n    title = {Randomized algorithms for the low-rank approximation of matrices},\n    journal = {Proceedings of the National Academy of Sciences},\n    volume = {104},\n    issue = {51},\n    year = {2007},\n    pages = {20167--20172},\n    doi  = {10.1073/pnas.0709640104}\n}\n\n\n\n"
},

{
    "location": "library/public.html#Norm-estimate-1",
    "page": "Public",
    "title": "Norm estimate",
    "category": "section",
    "text": "rnorm\nrnorms"
},

{
    "location": "library/public.html#IterativeSolvers.reig",
    "page": "Public",
    "title": "IterativeSolvers.reig",
    "category": "Function",
    "text": "reig(A, l)\n\nCompute the spectral (Eigen) decomposition of A using a randomized algorithm.\n\nArguments\n\nA: input matrix.\n\nl::Int: number of eigenpairs to find.\n\nOutput\n\n::Base.LinAlg.Eigen: eigen decomposition.\n\nImplementation note\n\nThis is a wrapper around eigfact_onepass() which uses the randomized samples found using srft(l).\n\nReferences\n\n@article{Halko2011,\n    author = {Halko, N and Martinsson, P G and Tropp, J A},\n    doi = {10.1137/090771806},\n    journal = {SIAM Review},\n    month = jan,\n    number = {2},\n    pages = {217--288},\n    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},\n    volume = {53},\n    year = {2011}\n}\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.rsvdfact",
    "page": "Public",
    "title": "IterativeSolvers.rsvdfact",
    "category": "Function",
    "text": "rsvdfact(A, n, p=0)\n\nCompute partial singular value decomposition of A using a randomized algorithm.\n\nArguments\n\nA: input matrix.\n\nn::Int: number of singular value/vector pairs to find.\n\np::Int=0: number of extra vectors to include in computation.\n\nOutput\n\n::SVD: singular value decomposition.\n\nWarning\n\nThis variant of the randomized singular value decomposition is the most commonly found implementation but is not recommended for accurate computations, as it often has trouble finding the n largest singular pairs, but rather finds n large singular pairs which may not necessarily be the largest.\n\nImplementation note\n\nThis function calls rrange, which uses naive randomized rangefinding to compute a basis for a subspace of dimension n (Algorithm 4.1 of cite{Halko2011}), followed by svdfact_restricted(), which computes the exact SVD factorization on the restriction of A to this randomly selected subspace (Algorithm 5.1 of cite{Halko2011}).\n\nAlternatively, you can mix and match your own randomized algorithm using any of the randomized range finding algorithms to find a suitable subspace and feeding the result to one of the routines that computes the SVD restricted to that subspace.\n\nReferences\n\n@article{Halko2011,\n    author = {Halko, N and Martinsson, P G and Tropp, J A},\n    doi = {10.1137/090771806},\n    journal = {SIAM Review},\n    month = jan,\n    number = {2},\n    pages = {217--288},\n    title = {Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions},\n    volume = {53},\n    year = {2011}\n}\n\n\n\n"
},

{
    "location": "library/public.html#IterativeSolvers.rsvd_fnkz",
    "page": "Public",
    "title": "IterativeSolvers.rsvd_fnkz",
    "category": "Function",
    "text": "rsvd_fnkz(A, k)\n\nCompute the randomized SVD by iterative refinement from randomly selected columns/rows.\n\nArguments\n\nA: matrix whose SVD is desired.\n\nk::Int: desired rank of approximation (k ≤ min(m, n)).\n\nKeywords\n\nl::Int = k: number of columns/rows to sample at each iteration (1 ≤ l ≤ k).\n\nN::Int = minimum(size(A)): maximum number of iterations.\n\nϵ::Real = prod(size(A))*eps(): relative threshold for convergence, as measured by growth of the spectral norm.\n\nmethod::Symbol = :eig: problem to solve.\n\n:eig: eigenproblem.\n:svd: singular problem.\n\nverbose::Bool = false: print convergence information at each iteration.\n\nOutput\n\nSVD object of rank ≤ k.\n\nReferences\n\n@inproceedings{,\n    author={Friedland, S. and Niknejad, A. and Kaveh, Mostafa and Zare, H.},\n    booktitle={System of Systems Engineering, 2006 IEEE/SMC International Conference on},\n    title={Fast Monte-Carlo low rank approximations for matrices},\n    year={2006},\n    month={April},\n    pages={218--223},\n    doi={10.1109/SYSOSE.2006.1652299}\n}\n\n\n\n"
},

{
    "location": "library/public.html#Randomized-singular-value-decomposition-1",
    "page": "Public",
    "title": "Randomized singular value decomposition",
    "category": "section",
    "text": "reig\nrsvdfact\nrsvd_fnkz"
},

{
    "location": "library/internal.html#",
    "page": "Internal",
    "title": "Internal",
    "category": "page",
    "text": ""
},

{
    "location": "library/internal.html#Internal-Documentation-1",
    "page": "Internal",
    "title": "Internal Documentation",
    "category": "section",
    "text": "Documentation for IterativeSolvers.jl's internals.Pages = [\"internal.md\"]\nDepth = 4"
},

{
    "location": "library/internal.html#Index-1",
    "page": "Internal",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internal.md\"]"
},

{
    "location": "library/internal.html#ConvergenceHistory-Internals-1",
    "page": "Internal",
    "title": "ConvergenceHistory Internals",
    "category": "section",
    "text": "TypealiasesIterativeSolvers.PlainHistory\nIterativeSolvers.RestartedHistoryFunctionsIterativeSolvers.nextiter!\nIterativeSolvers.reserve!\nIterativeSolvers.shrink!\nIterativeSolvers.setmvps\nIterativeSolvers.setmtvps\nIterativeSolvers.setconv\nIterativeSolvers.showplot"
},

{
    "location": "library/internal.html#IterativeSolvers.lastvec",
    "page": "Internal",
    "title": "IterativeSolvers.lastvec",
    "category": "Function",
    "text": "lastvec(K)\n\nGet last vector computed in the Krylov subspace K.\n\nArguments\n\nK::KrylovSubspace: Krylov subspace.\n\nOutput\n\n::Vector: last vector.\n\n\n\n"
},

{
    "location": "library/internal.html#IterativeSolvers.nextvec",
    "page": "Internal",
    "title": "IterativeSolvers.nextvec",
    "category": "Function",
    "text": "nextvec(K)\n\nCompute next vector in the Krylov subspace K.\n\nArguments\n\nK::KrylovSubspace: Krylov subspace.\n\nOutput\n\n::Vector: next vector.\n\n\n\n"
},

{
    "location": "library/internal.html#IterativeSolvers.init!",
    "page": "Internal",
    "title": "IterativeSolvers.init!",
    "category": "Function",
    "text": "init!(K, v)\n\nInitialize the KrylovSubspace K with a specified nonunit vector.\n\nArguments\n\nK::KrylovSubspace: Krylov subspace.\n\nv::Vector: initial vector.\n\n\n\n"
},

{
    "location": "library/internal.html#IterativeSolvers.initrand!-Tuple{IterativeSolvers.KrylovSubspace}",
    "page": "Internal",
    "title": "IterativeSolvers.initrand!",
    "category": "Method",
    "text": "initrand!(K)\n\nInitialize the Krylov subspace K with a random unit vector.\n\nArguments\n\nK::KrylovSubspace: Krylov subspace.\n\n\n\n"
},

{
    "location": "library/internal.html#IterativeSolvers.appendunit!",
    "page": "Internal",
    "title": "IterativeSolvers.appendunit!",
    "category": "Function",
    "text": "appendunit!(K, w)\n\nAppend normalize(w) vector to krylob subspace K.\n\nArguments\n\nK::KrylovSubspace: Krylov subspace.\n\nw::Vector: vector to append.\n\n\n\n"
},

{
    "location": "library/internal.html#IterativeSolvers.orthogonalize",
    "page": "Internal",
    "title": "IterativeSolvers.orthogonalize",
    "category": "Function",
    "text": "orthogonalize{T}(v, K, p)\n\nOrthogonalize a vector v against the last p basis vectors defined by the Krylov subspace K.\n\nArguments\n\nv::Vector: vector to orthogonalize.\n\nK::KrylovSubspace: subspace to orthogonalize v against.\n\np::Int=K.order: last p remembered vectors.\n\nKeywords\n\nmethod::Symbol=:ModifiedGramSchmidt: orthogonalization method. Choose from:\n\n:GramSchmidt: Gram Schmidt method.\n:ModifiedGramSchmidt: Modified Gram Schmidt method.\n:Householder: Householder method.\n\nnormalize::Bool=false: normalize vector.\n\nOutput\n\n::Vector: orthogonalized vector.\n\n::Vector: orthogonalization coefficients.\n\n\n\n"
},

{
    "location": "library/internal.html#KrylovSubspace-Internals-1",
    "page": "Internal",
    "title": "KrylovSubspace Internals",
    "category": "section",
    "text": "FunctionsIterativeSolvers.lastvec\nIterativeSolvers.nextvec\nIterativeSolvers.init!\nIterativeSolvers.initrand!(::KrylovSubspace)\nIterativeSolvers.appendunit!\nIterativeSolvers.orthogonalize"
},

{
    "location": "library/internal.html#Other-Functions-1",
    "page": "Internal",
    "title": "Other Functions",
    "category": "section",
    "text": "IterativeSolvers.idfact\nIterativeSolvers.isconverged\nIterativeSolvers.thickrestart!\nIterativeSolvers.harmonicrestart!\nIterativeSolvers.plot_collection\nIterativeSolvers.plotable\nIterativeSolvers.Adivtype\nIterativeSolvers.Amultype\nIterativeSolvers.randx\nIterativeSolvers.zerox\nIterativeSolvers.update!\nIterativeSolvers.initrand!(::Vector)"
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
    "text": "Contributions are always welcome, as are feature requests and suggestions. Feel free to open issues and pull requests at any time. If you aren't familiar with git or Github please start now.It is important to note that almost every method in the package has documentation, to know what it does simply use ?<method> in the terminal.julia> using IterativeSolvers\n\nhelp?> IterativeSolvers.Adivtype\n  Adivtype(A, b)\n\n  Determine type of the division of an element of b against an element of A:\n\n  typeof(one(eltype(b))/one(eltype(A)))"
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
    "text": "Each iterative method method must log information using the inner ConvergenceHistory type. When information is not necessary to be stored (plot is set to false) then instead of ConvergenceHistory create a DummyHistory, this type has the same calls ConvergenceHistory does but without storing anything.There are two types of ConvergenceHistory: plain and restarted. The only real difference between the two is how they are plotted and how the number of restarts is calculated, everything else is the same.Before logging information space must always be reserved.log = ConvergenceHistory()\nlog[:tol] = tol\nreserve!(log,:betas, maxiter) # Vector of length maxiter\nreserve!(log,:conv, maxiter, T=BitArray) # Vector of length maxiter\nreserve!(log,:ritz, maxiter, k) # Matrix of size (maxiter, k)To store information at each iteration use push!.push!(log, :conv, conv)\npush!(log, :ritz, F[:S][1:k])\npush!(log, :betas, L.β)To advance the log index to the next iteration use nextiter!.nextiter!(log)A more detailed explanation of all the functions is in both the public and internal documentation of ConvergenceHistory.The most rich example of the usage of ConvergenceHistory is in svdl."
},

{
    "location": "about/CONTRIBUTING.html#Adding-benchmarks-1",
    "page": "Contributing",
    "title": "Adding benchmarks",
    "category": "section",
    "text": "The Benchmarks tab of the documentation is built automatically with Travis. Any benchmark added will be displayed automatically after a successful pull request.The benchmark suite gets built doing a cross product between the available matrices and available methods, if there are n methods and m linear operators then n*m will be the upper limit of benchmarks to be made. Some methods are not compatible with certain matrices, to avoid generating unnecessary benchmarks each method and matrix has traits, linear operator traits are inspired from MatrixDepot.jl.Method traitsaccessible : Method accesses the linear operator's fields.\ninverse    : A's Inverse must exist.\nsymmetric  : A's must be symmetric.\npos-def    : A's must be definite.Linear Operator traitsaccessible : Is accessible.\ninverse    : A is exist.\nsymmetric  : A is symmetric.\npos-def    : A is definite.\neigen      : Part of the eigensystem of the matrix is explicitly known.\ngraph      : An adjacency matrix of a graph.\nill-cond   : The matrix is ill-conditioned for some parameter values.\nrandom     : The matrix has random entries.\nregprob    : The output is a test problem for Regularization Methods.\nsparse     : The matrix is sparse.A benchmark between a method and a linear operator will be made if and only if the traits of the method is subset of the traits of the linear operator.Benchmarks are stored in Benchmarks.jl. To add a method use addEqMethod.addEqMethod(methods, \"jacobi\", jacobi, [\"inverse\",\"accessible\"])\naddEqMethod(methods, \"gauss_seidel\", gauss_seidel, [\"inverse\",\"accessible\"])\naddEqMethod(methods, \"sor\", sor, [\"inverse\",\"accessible\"])\naddEqMethod(methods, \"ssor\", ssor, [\"inverse\",\"accessible\", \"symmetric\"])\naddEqMethod(methods, \"cg\", cg, [\"inverse\", \"symmetric\", \"pos-def\"])\naddEqMethod(methods, \"gmres\", gmres, [\"inverse\"])\naddEqMethod(methods, \"lsqr\", lsqr, [\"inverse\"])\naddEqMethod(methods, \"chebyshev\", chebyshev, [\"inverse\", \"accessible\"])Here methods is a dictionary, the second argument is the name to be displayed in the benchmarks, the third argument is the function and the fourth is the traits. Every function has a predetermined call in buildCall function.To add an equation use addEquation.#Sparse matrix equations\naddEquation(\n    equations, \"Poisson\", [\"Sparse\", \"Poisson\"],\n    [\"sparse\",\"inverse\", \"symmetric\", \"pos-def\", \"eigen\", \"accessible\"],\n    :(matrixdepot(\"poisson\",4))\n    )\n\n#Function matrix equations\naddEquation(\n    equations, \"SOLtest\", [\"Function\", \"SOLtest\"],\n    [\"function\",\"inverse\"],\n    :(buildSol(10)),\n    10\n)Here equations is a dictionary, the second argument is the name to be displayed in the benchmarks, the third argument is the path inside the BenchmarkGroup type the fourth argument is the traits, the fifth is the matrix generator and the sixth is the size of the matrix. The size of the matrix has to be passed when it is impossible to deduce the dimension from the generator, in this case buildSol generates a function and not a matrix.To add a custom benchmark use directly the suite variable which is the BenchmarkGroup of the package, to know more of this type check BenchmarkTools.jl."
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
