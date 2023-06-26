"""
    BipartiteProjectionFailed

This exception is thrown when a unipartite network cannot be projected into a
bipartite network. This happens when some species have both successors *and*
predecessors.
"""
struct BipartiteProjectionFailed <: Exception
end

function Base.showerror(io::IO, e::BipartiteProjectionFailed)
    message = "The network cannot be made unipartite as some of its species overlap"
    return print(io, message)
end