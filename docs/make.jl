using Documenter, IterativeSolvers

makedocs(
	modules = [IterativeSolvers],
	format = :html,
	sitename = "IterativeSolvers.jl",
	pages = [
		"Home" => "index.md",
		"Getting started" => "getting_started.md",
		"Preconditioning" => "preconditioning.md",
		"Linear systems" => [
			"Conjugate Gradients" => "linear_systems/cg.md",
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
		],
		"SVDL" => "svd/svdl.md",
		"Randomized algorithms" => "randomized.md",
		"The iterator approach" => "iterators.md",
		# "Additional resources" => [
		# 	# "Public" => "library/public.md",
		# 	# "Internal" => "library/internal.md",
		# ],
		"About" => [
			"Contributing" => "about/CONTRIBUTING.md",
			"License" => "about/license.md",
	    ]
	]
)

deploydocs(
	repo = "github.com/JuliaMath/IterativeSolvers.jl.git",
	target = "build",
	osname = "linux",
	julia  = "0.6",
	deps = nothing,
	make = nothing,
)
