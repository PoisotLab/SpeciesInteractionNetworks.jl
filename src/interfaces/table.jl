Tables.istable(::Type{<:SpeciesInteractionNetwork}) = true
Tables.rowaccess(::Type{<:SpeciesInteractionNetwork}) = true
Tables.rows(N::SpeciesInteractionNetwork) = interactions(N)
Tables.schema(::SpeciesInteractionNetwork{<:Partiteness{T1}, <:Interactions{T2}}) where {T1, T2} = Tables.Schema((:from, :to, :value), (T1, T1, T2))