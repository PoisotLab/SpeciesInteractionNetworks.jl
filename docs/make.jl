push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterVitepress
using DocumenterCitations
using SpeciesInteractionNetworks
using Literate

bibliography = CitationBibliography(joinpath(@__DIR__, "SpeciesInteractionNetworks.bib"), style=:authoryear)

doc_src = first(splitdir(@__FILE__))

usecases = readdir(joinpath(doc_src, "src/use_cases/"); join=true)
literate_config = Dict("credit" => false)
filter!(endswith(".jl"), usecases)
for usecase in usecases
    Literate.markdown(usecase, first(splitdir(usecase)); config=literate_config)
end

makedocs(
    bibliography;
    sitename="SpeciesInteractionNetworks",
    authors="Timothée Poisot",
    modules=[SpeciesInteractionNetworks],
    format=MarkdownVitpress(
        repo="github.com/PoisotLab/SpeciesInteractionNetworks.jl",
    ),
)

deploydocs(;
    repo="github.com/PoisotLab/SpeciesInteractionNetworks.jl.git",
    devbranch="main",
    push_preview=true,
)
