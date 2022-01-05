using Documenter, IterativeSolvers

makedocs(
	modules = [IterativeSolvers],
	format = Documenter.HTML(
		# Disable pretty URLs during manual testing
		prettyurls = get(ENV, "CI", nothing) == "true",
		# Set canonical URL to GitHub pages URL
		canonical = "https://julialinearalgebra.github.io/IterativeSolvers.jl/stable"
  ),
	doctest = false,
	clean = true,
	sitename = "IterativeSolvers.jl",
	pages = [
		"Home" => "index.md",
		"Getting started" => "getting_started.md",
		"Preconditioning" => "preconditioning.md",
		"Linear systems" => [
			"Conjugate Gradient" => "linear_systems/cg.md",
			"Chebyshev iteration" => "linear_systems/chebyshev.md",
			"MINRES" => "linear_systems/minres.md",
			"BiCGStab(l)" => "linear_systems/bicgstabl.md",
			"IDR(s)" => "linear_systems/idrs.md",
			"Restarted GMRES" => "linear_systems/gmres.md",
			"LSMR" => "linear_systems/lsmr.md",
			"LSQR" => "linear_systems/lsqr.md",
			"Stationary methods" => "linear_systems/stationary.md"
		],
		"Eigenproblems" => [
			"Power method" => "eigenproblems/power_method.md",
			"LOBPCG" => "eigenproblems/lobpcg.md"
		],
		"SVDL" => "svd/svdl.md",
		"The iterator approach" => "iterators.md",
		"About" => [
			"Contributing" => "about/CONTRIBUTING.md",
			"License" => "about/license.md",
	    ]
	]
)

deploydocs(repo = "github.com/JuliaLinearAlgebra/IterativeSolvers.jl.git")
