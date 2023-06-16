"""
    StructuralModel

...
"""
abstract type StructuralModel end

"""
    NicheModel

The niche model of food webs is one of the most emblematic food web models.
Based on a given species richness and *expected* connectance, it generates food
webs with properties that closely reproduce empirical networks. One particular
property of the food webs generated with this model is that they are *interval*,
which is to say that there are not gaps in the diet of species.

###### References

[Williams2000Simple](@citet*)
"""
abstract type NicheModel <: StructuralModel end
