function _slurp(network::Mangal.MangalNetwork, query::Pair...)
    page_length = 50
    result_size = count(Mangal.MangalInteraction, network, query...)
    network_content = Vector{Mangal.MangalInteraction}(undef, result_size)
    total_pages = ceil(Int64, result_size / page_length)
    for current_page in 1:total_pages
        this_page_query = ["page" => "$(current_page - 1)", "count" => page_length, query...]
        offset_start = (current_page - 1) * page_length + 1
        offset_end = min(offset_start + page_length - 1, result_size)
        network_content[offset_start:offset_end] .= Mangal.interactions(network, this_page_query...)
    end
    return network_content
end

"""
    mangalnetwork(MN::Mangal.MangalNetwork, query::Pair...; taxonlevel::Bool=false)

Retrieves a network from [mangal.io](https://mangal.io), using the `Mangal.jl`
API wrapper. The keyword argument `taxonlevel` specifies whether the original
*node* or its *reference taxon* must be used. Note that not all nodes have an
unambiguously attached reference taxon.

###### References

[Poisot2016mangal](@citet*)
"""
function mangalnetwork(MN::Mangal.MangalNetwork, query::Pair...; taxonlevel::Bool=false)
    edgelist = _slurp(MN, query...)
    nodelist = [e.from for e in edgelist] ∪ [e.to for e in edgelist]
    if taxonlevel
        nodelist = unique(filter(!ismissing, map(n -> n.taxon, nodelist)))
    end
    nodes = Unipartite(nodelist)
    T = typeof(edgelist[1].strength)
    edges = Quantitative(zeros(T, size(nodes)))
    N = SpeciesInteractionNetwork(nodes, edges)
    for edge in edgelist
        rf = taxonlevel ? edge.from.taxon : edge.from
        rt = taxonlevel ? edge.to.taxon : edge.to
        rq = edge.strength
        if !(ismissing(rf)|ismissing(rt))
            N[rf, rt] = rq
        end
    end
    return N
end