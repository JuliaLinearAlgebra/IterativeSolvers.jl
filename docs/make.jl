using Documenter, IterativeSolvers

makedocs(
	modules = [IterativeSolvers],
	format = :html,
	sitename = "IterativeSolvers.jl",
	pages = [
		"Home" => "index.md",
		"Manual" => "user_manual.md",
		"Linear systems" => [
			"Conjugate Gradients" => "library/cg.md",
			"Chebyshev iteration" => "library/chebyshev.md",
			"MINRES" => "library/minres.md",
			"BiCGStab(l)" => "library/bicgstabl.md",
			"Restarted GMRES" => "library/gmres.md",
			"Stationary methods" => "library/stationary.md"
		],
		"Eigenproblems" => [
			"Power method" => "library/power_method.md",
		],
		"Preconditioning" => "preconditioning.md",
		"Library" => [
			"Public" => "library/public.md",
			"Internal" => "library/internal.md",
		],
		"About" => [
		"Contributing" => "about/CONTRIBUTING.md",
		"License" => "about/license.md",
	    ]
	]
)

deploydocs(
	repo = "github.com/JuliaMath/IterativeSolvers.jl",
	target = "build",
	deps = nothing,
	make = nothing,
)
