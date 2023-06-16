# Building a network

!!! abstract

    Before getting started with the package itself, we will see how we can build a network, access its content, and iterate over the interactions. This page is intended to give you some intuitions about how the type system works, before reading more of the manual.

To begin with, we will load the package:

```@example 1
using SpeciesInteractionNetworks
```

## Creating a list of species

We will create a very small network, made of four species and their
interactions. The first step is to define a list of species:

```@example 1
species = [:fox, :vole, :hawk, :turnip]
```

In order to make sure that we are explicit about the type of network we are
working with, we will create a representation of this list of species that is
unipartite, using the [`Unipartite`](@ref) constructor:

```@example 1
nodes = Unipartite(species)
```

Note that the package is not considering information about the *ecological
nature* of the interaction, only (i) the structure of the community, as captured
by its [`Partiteness`](@ref), and later on about its [`Interactions`](@ref).

## Creating interactions

As with species, we want to represent interactions in a way that captures
ecological information. In this case, we will use binary interactions (true/0),
and work from a matrix, where the rows are the source of the interaction, and
the column is its destination. It means that interactions go *from* predator
*to* preys.

```@example 1
int_matrix = Bool[
    0 1 0 0;
    0 0 0 1;
    0 1 0 0;
    0 0 0 0
]
```

!!! info "About interaction as matrices"
    
    By specifying interactions as a matrix, it is fundamental that columns and orders are in the correct order. There are alternative ways to specify networks that do not rely on matrices (using tuples or pairs), but because most species interaction network data are represented as matrices, this is supported by the package.

As this network is binary, we will wrap this matrix into a [`Binary`](@ref)
collection of interactions:

```@example 1
edges = Binary(int_matrix)
```

## Assembling the network

The network itself is a collection of nodes and edges. There are a number of
specific checks performed when creating the network, to ensure that we cannot
create an object that makes no sense. These checks are done when calling
[`SpeciesInteractionNetwork`](@ref), which is the main type around which the
package is built.

```@example 1
network = SpeciesInteractionNetwork(nodes, edges)
```

The networks are iterable, *i.e.* we can walk through them, specifically one
interaction at a time:

```@example 1
for interaction in network
    println(interaction)
end
```

Internally, this is done by only returning the pairs of species that do not have
a value of zero. There is a way to capture all of the interactions at a time,
using [`interactions`](@ref):

``` @example 1
interactions(network)
```

## Basics of network exploration

We can also get a list of the species that establish an interaction with a given
species (in this case, predators):

```@example 1
predecessors(network, :vole)
```

Or the species witch which a given species establishes interactions (in this
case, preys):

```@example 1
successors(network, :fox)
```

Furthermore, we can return a subset (or more accurately a [`subgraph`](@ref)) of
the network, by giving a list of species:

```@example 1
interactions(subgraph(network, [:fox, :vole, :turnip]))
```

## Networks are editable

The content of networks can be modified. For example, to circumvent the issue of
needing to write the interaction matrix in the correct order, we can start with
an empty network:

```@example 1
netsize = (richness(nodes,1), richness(nodes, 2))
edges2 = Binary(zeros(Bool, netsize))
network2 = SpeciesInteractionNetwork(nodes, edges2)
interactions(network2)
```

We can then add the interactions one by one:

```@example 1
for interaction in [(:fox, :vole), (:hawk, :vole), (:vole, :turnip)]
    network2[interaction...] = true
end
interactions(network2)
```

## Networks are tables

All of the networks can be converted to a tabular data, for use with *e.g.* the
`DataFrames` package:

```@example 1
import DataFrames
DataFrames.DataFrame(network)
```
