using Documenter, IterativeSolvers

makedocs(modules = [IterativeSolvers])

deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-material"),
    repo = "github.com/lopezm94/IterativeSolvers.jl"
)
