push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterMarkdown
using QuantumCitations
using SpeciesInteractionNetworks

bibliography = CitationBibliography(joinpath(@__DIR__, "SpeciesInteractionNetworks.bib"), style = :authoryear)

makedocs(
    bibliography;
    sitename = "SpeciesInteractionNetworks",
    authors = "Timothée Poisot",
    modules = [SpeciesInteractionNetworks],
    format = Markdown(),
)

deploydocs(;
    deps = Deps.pip("mkdocs", "pygments", "python-markdown-math", "mkdocs-material", "mkdocs-git-revision-date-localized-plugin"),
    repo = "github.com/PoisotLab/SpeciesInteractionNetworks.jl.git",
    devbranch = "main",
    make = () -> run(`mkdocs build`),
    target = "site",
    push_preview = true,
)
