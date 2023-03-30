using EllipseSampling
using Documenter

DocMeta.setdocmeta!(EllipseSampling, :DocTestSetup, :(using EllipseSampling); recursive=true)

makedocs(;
    modules=[EllipseSampling],
    authors="JoelTrent <79883375+JoelTrent@users.noreply.github.com> and contributors",
    repo="https://github.com/JoelTrent/EllipseSampling.jl/blob/{commit}{path}#{line}",
    sitename="EllipseSampling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JoelTrent.github.io/EllipseSampling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JoelTrent/EllipseSampling.jl",
    devbranch="main",
)
