module SpeciesInteractionNetworks

# Dependencies
using Combinatorics
using DataStructures
using DelimitedFiles
using Distributions
using LinearAlgebra
using Random
using Requires
using SparseArrays
using Statistics
using StatsBase
using TestItems
import Tables
import Mangal
import Graphs

# Various utilities for probabilities
# include(joinpath(".", "misc", "probabilities.jl"))

# Tests to define what can be used in base types
# include(joinpath(".", "misc", "init_tests.jl"))

include("types/exceptions.jl")
export BipartiteProjectionFailed

include("types/declarations.jl")
export Partiteness, Interactions
export Bipartite, Unipartite
export Binary, Quantitative, Probabilistic
export SpeciesInteractionNetwork

include("types/constructors.jl")

include("types/interface.jl")
export species, richness

include("types/copy.jl")

include("types/render.jl")
export render

include("interfaces/abstractarray.jl")
include("interfaces/iteration.jl")
include("interfaces/table.jl")
include("interfaces/broadcast.jl")
include("interfaces/linearalgebra.jl")
include("interfaces/graphs.jl")

export svd, rank, diag
export complexity, tsvd, rdpg

include("basics/interactions.jl")
export interactions

include("basics/neighbors.jl")
export successors, predecessors, neighbors

include("basics/subgraphs.jl")
export subgraph

include("random/permutations.jl")
export PermutationConstraint, Degree, Generality, Vulnerability, Connectance
export swap!

include("random/structural.jl")
export StructuralModel

include("random/structural/nichemodel.jl")
export NicheModel

export structuralmodel

include("random/draws.jl")
export randomdraws

include("random/linearfilter.jl")
export linearfilter, nullmodel

include("random/contributions.jl")
export speciescontribution

include("micro_level/degree.jl")
export degree
export generality, vulnerability

include("micro_level/specificity.jl")
export specificity

include("micro_level/centrality.jl")
export CentralityMethod
export KatzCentrality, EigenvectorCentrality
export ClosenessCentrality, ResidualClosenessCentrality, GeneralizedClosenessCentrality
export centrality

include("meso_level/paths.jl")
export ShortestPathMethod
export normalize
export shortestpath, pathbetween

include("micro_level/foodwebs.jl")
export distancetobase

include("meso_level/paths/BellmanFord.jl")
export BellmanFord

include("meso_level/paths/Dijkstra.jl")
export Dijkstra

include("meso_level/motifs/generators.jl")
export motifs

include("meso_level/motifs/permutations.jl")
include("meso_level/motifs.jl")
export findmotif

include("meso_level/modularity.jl")
export modularity

include("meso_level/modularity/labelpropagation.jl")
export labelpropagation

include("macro_level/connectance.jl")
export connectance, links, linkagedensity
export connectance_variance, links_variance, linkagedensity_variance

include("macro_level/eta.jl")
export η

include("macro_level/nodf.jl")
export nodf

include("macro_level/spectralradius.jl")
export spectralradius

include("macro_level/degeneracy.jl")
export isdegenerate, isdisconnected
export simplify

include("meta_level/set.jl")
include("meta_level/partitions.jl")
export BetaDivComponent
export βS, βOS, βWN
export betadiversity

include("meta_level/measures.jl")
export KGL01, KGL02, KGL03, KGL04, KGL05, KGL06, KGL07, KGL08, KGL09, KGL10, KGL11, KGL12, KGL13, KGL14, KGL15, KGL16, KGL17, KGL18, KGL19, KGL20, KGL21, KGL22, KGL23, KGL24

include("data/mangal.jl")
export mangalnetwork

# include(joinpath(".", "types", "conversions.jl"))

# Datasets
# include(joinpath(".", "misc", "data.jl"))
# export web_of_life, nz_stream_foodweb, pajek

#=
function __init__()
   @require Mangal="b8b640a6-63d9-51e6-b784-5033db27bef2" begin
      _check_species_validity(::Mangal.MangalReferenceTaxon) = nothing
      _check_species_validity(::Mangal.MangalNode) = nothing
   end
   @require GBIF="ee291a33-5a6c-5552-a3c8-0f29a1181037" begin
      _check_species_validity(::GBIF.GBIFTaxon) = nothing
   end
   @require NCBITaxonomy="f88b31d2-eb98-4433-b52d-2dd32bc6efce" begin
      _check_species_validity(::NCBITaxonomy.NCBITaxon) = nothing
   end
end
=#

# General useful manipulations
# include(joinpath(".", "utilities", "comparisons.jl"))
# include(joinpath(".", "utilities", "overloads.jl"))
# include(joinpath(".", "utilities", "utilities.jl"))
# export species, interactions, has_interaction, richness, nodiagonal, nodiagonal!, adjacency

# Degree
# include(joinpath(".", "links", "degree.jl"))
# export degree, degree_var

# include(joinpath(".", "links", "connectance.jl"))
# export links, L, links_var, connectance, connectance_var,
#    linkage_density, link_number

# Random networks from structural models
# include(joinpath(".", "structuralmodels", "cascademodel.jl"))
# export cascademodel
# include(joinpath(".", "structuralmodels", "mpnmodel.jl"))
# export mpnmodel
# include(joinpath(".", "structuralmodels", "nestedhierarchymodel.jl"))
# export nestedhierarchymodel

# Centrality
# include(joinpath(".", "community", "centrality.jl"))
# export centrality_katz, centrality_degree, centrality_eigenvector
# export centrality_closeness, centrality_harmonic

# Motifs
# include(joinpath(".", "community", "motifs.jl"))
# export find_motif, expected_motif_count, unipartitemotifs

# Overlap
# include(joinpath(".", "community", "overlap.jl"))
# export overlap
# export AJS, EAJS

# Overlap
# include(joinpath(".", "community", "resilience.jl"))
# export resilience
# export symmetry, heterogeneity
# export s
# export σ_in, σ_out

# Modularity
# include(joinpath(".", "modularity", "utilities.jl"))
# export Q, Qr
# 
# include(joinpath(".", "modularity", "labelpropagation.jl"))
# export lp, salp
# 
# include(joinpath(".", "modularity", "starters.jl"))
# export n_random_modules, each_species_its_module
# 
# include(joinpath(".", "modularity", "brim.jl"))
# export brim
# 
# include(joinpath(".", "modularity", "roles.jl"))
# export functional_cartography

# Food webs
# include(joinpath(".", "foodwebs", "trophiclevels.jl"))
# export distance_to_producer, trophic_level, omnivory

# include(joinpath(".", "information", "entropy.jl"))
# export entropy, make_joint_distribution, mutual_information, conditional_entropy,
#    variation_information, diff_entropy_uniform, information_decomposition,
#    convert2effective, potential_information

# include(joinpath(".", "information", "otsin.jl"))
# export optimaltransportation

end
