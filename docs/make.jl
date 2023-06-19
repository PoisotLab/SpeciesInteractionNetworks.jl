push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterMarkdown
using QuantumCitations
using SpeciesInteractionNetworks
using Literate

bibliography = CitationBibliography(joinpath(@__DIR__, "SpeciesInteractionNetworks.bib"), style = :authoryear)

usecases = readdir("src/use_cases/"; join=true)
literate_config = Dict(credit => false)
filter!(endswith(".jl"), usecases)
for usecase in usecases
    Literate.markdown(usecase, outputdir=first(splitdir(usecase)); config=literate_config)
end

makedocs(
    bibliography;
    sitename = "SpeciesInteractionNetworks",
    authors = "Timothée Poisot",
    modules = [SpeciesInteractionNetworks],
    format = Markdown(),
)

deploydocs(;
    deps = Deps.pip("mkdocs", "pygments", "python-markdown-math", "mkdocs-material"),
    repo = "github.com/PoisotLab/SpeciesInteractionNetworks.jl.git",
    devbranch = "main",
    make = () -> run(`mkdocs build`),
    target = "site",
    push_preview = true,
)
