function pdi(impacts::AbstractArray{T,1}) where {T <: Number}
    norm_impact = impacts/maximum(impacts)
    paired_differences = (1.0 .- sort(norm_impact/maximum(norm_impact), rev=true))[2:end]
    return sum(paired_differences)/length(paired_differences)
end

@testitem "The PDI of a perfect specialist is one" begin
    @test SpeciesInteractionNetworks.pdi([true, false, false]) == 1.0
    @test SpeciesInteractionNetworks.pdi([2, 0, 0]) == 1.0
    @test SpeciesInteractionNetworks.pdi([1.0, 0.0, 0.0]) == 1.0
end

@testitem "The PDI of a perfect generalist is zero" begin
    @test SpeciesInteractionNetworks.pdi([true, true, true]) == 0.0
    @test SpeciesInteractionNetworks.pdi([2, 2, 2]) == 0.0
    @test SpeciesInteractionNetworks.pdi([1.0, 1.0, 1.0]) == 0.0
end

"""
    specificity(N::SpeciesInteractionNetwork{<:Partiteness, <:Union{Binary,Quantitative}})

For a deterministic network, this function will return a dictionary mapping each
top-level species to its specificity. The same index (Paired Differences Index)
is used for binary and quantitative networks. A value of one corresponds to
maximum specialism, and a value of zero to maximum generalism. The index is
symmetrical, so that a species with a value of one half is neither specialist
nor generalist.

###### References

[Poisot2012comparative](@cite)
"""
function specificity(N::SpeciesInteractionNetwork{<:Partiteness, <:Union{Binary,Quantitative}})
    return Dict([s => pdi(N[s,:]) for s in species(N, 1)])
end