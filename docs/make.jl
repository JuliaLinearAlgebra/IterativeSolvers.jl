using Documenter, IterativeSolvers

#include("buildmarkdown.jl")

makedocs(modules = [IterativeSolvers])

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math", "mkdocs-material"),
    repo = "github.com/lopezm94/IterativeSolvers.jl"
)
