using ExtremalOptimization
using Documenter

makedocs(;
    modules=[ExtremalOptimization],
    authors="Francesco Alemanno <francescoalemanno710@gmail.com> and contributors",
    repo="https://github.com/francescoalemanno/ExtremalOptimization.jl/blob/{commit}{path}#L{line}",
    sitename="ExtremalOptimization.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://francescoalemanno.github.io/ExtremalOptimization.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/francescoalemanno/ExtremalOptimization.jl",
)
