using Revise
using EllipseSampling
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(EllipseSampling, :DocTestSetup, :(using EllipseSampling); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

makedocs(;
    modules=[EllipseSampling],
    authors="JoelTrent <79883375+JoelTrent@users.noreply.github.com> and contributors",
    sitename="EllipseSampling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JoelTrent.github.io/EllipseSampling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quick Start" => "quick_start.md",
        "User Interface" => "user_interface.md",
        "Internal Library" => "internal_library.md",
        "References" => "references.md"
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/JoelTrent/EllipseSampling.jl",
    devbranch="main",
)
