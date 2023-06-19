import Pkg

# Switch the project to the documentation project
Pkg.activate("docs/")

# Add the main package (dev)
Pkg.develop(Pkg.PackageSpec(;path="."))

# Instantiate everything
Pkg.instantiate()

# Build the docs
include("docs/make.jl")

# Drop into the docs folder and run mkdocs
cd("docs")
run(`mkdocs build`)

# Drop back to the main folder and remove the dev package
cd("..")
Pkg.rm("SpeciesInteractionNetworks")