using Documenter
using Radix

DocMeta.setdocmeta!(Radix, :DocTestSetup, :(using Radix); recursive=true)

makedocs(;
    modules=[Radix],
    authors="Paul Barrett <pebarrett@gmail.com>",
    sitename="Radix.jl",
    format=Documenter.HTML(;
        canonical="https://juliaastro.org/Radix.jl",
        edit_link="main",
        assets=String[],
        size_threshold_warn=500 * 2^10,
        size_threshold=2^20,
    ),
    pages=[
        "Home" => "index.md",
        "Using Radix" => "usage.md",
        "Reference" => "api.md",
        "Field notes" => "field-names.md",
    ],
    warnonly=[:missing_docs, :cross_references, :autodocs_block],
)

deploydocs(;
    repo="github.com/JuliaAstro/Radix.jl",
    devbranch="main",
)
