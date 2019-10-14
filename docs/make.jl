using Documenter, AnimatedOptimization

makedocs(;
    modules=[AnimatedOptimization],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/schrimpf/AnimatedOptimization.jl/blob/{commit}{path}#L{line}",
    sitename="AnimatedOptimization.jl",
    authors="Paul Schrimpf <paul.schrimpf@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/schrimpf/AnimatedOptimization.jl",
)
