using Documenter, IterativeSolvers

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-bootstrap"),
    repo = "github.com/lopezm94/IterativeSolvers.jl",
    julia  = "0.4.5"
)
