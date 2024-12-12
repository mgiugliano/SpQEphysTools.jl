using SpQEphysTools
using Documenter

DocMeta.setdocmeta!(SpQEphysTools, :DocTestSetup, :(using SpQEphysTools); recursive=true)

makedocs(;
    modules=[SpQLib],
    authors="Michele GIUGLIANO",
    sitename="SpQEphysTools.jl",
    format=Documenter.HTML(;
        canonical="https://mgiugliano.github.io/SpQEphysTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mgiugliano/SpQEphysTools.jl",
    devbranch="main",
)
