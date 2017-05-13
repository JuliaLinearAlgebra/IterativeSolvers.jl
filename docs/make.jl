using Documenter, IterativeSolvers

makedocs(
	modules = [IterativeSolvers],
	format = :html,
	sitename = "IterativeSolvers.jl",
	pages = [
		"Home" => "index.md",
		"Manual" => "user_manual.md",
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
