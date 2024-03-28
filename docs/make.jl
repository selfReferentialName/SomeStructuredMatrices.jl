using SomeStructuredMatrices
using Documenter

DocMeta.setdocmeta!(SomeStructuredMatrices, :DocTestSetup, :(using SomeStructuredMatrices); recursive=true)

makedocs(;
    modules=[SomeStructuredMatrices],
    authors="Jenifer De Jager (jraind@vt.edu) and contributors",
    sitename="SomeStructuredMatrices.jl",
    format=Documenter.HTML(;
        canonical="https://selfReferentialName.github.io/SomeStructuredMatrices.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/selfReferentialName/SomeStructuredMatrices.jl",
    devbranch="main",
)
