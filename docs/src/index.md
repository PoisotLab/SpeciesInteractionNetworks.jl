# Species Interaction Networks

The `SpeciesInteractionNetworks` package enables analyses of species interaction
networks in *Julia*. The list of implemented measures (and more broadly, the
design of the package) closely follows the recommendations in
[Delmas2019Analysing](@citet).

## Organisation of the package

The measures in the documentation are organized by level of organisation,
following a simple convention.

*Micro*-level measures return information at the level of the node. Typically,
these will return a value for each species in the network. When called on an
entire network, these will return a dictionary mapping the species to their
value; when called on a network *and* a species, these will return the value of
the species.

*Meso*-level measures return information about more than one species, but not
about the entire network. Examples are shortest paths, motif enumeration, etc.
The type of returned objects for these functions is more diverse and
domain-specific.

*Macro*-level measures are summaries of the entire network, and typically return
one value for the entire network.

Finally, *meta*-level measures are intended to be applied to a collection of
networks; network dissimilarity measures are one example.

## Why a new package?

TL;DR: technical bankruptcy.

The very first lines of
[`EcologicalNetworks`](https://github.com/PoisotLab/EcologicalNetworks.jl/) were
written for *Julia* 0.4 (or 0.5?), and the language has changed immensely since
then. By working with the package consistently over the years, it became
apparent that we were often working *around* design decisions. A lot of code was
also awkward to maintain, in part because it dated from before changes in the
langage.

It wasn't broken, but it wasn't *boss*. Think of this package as a re-telling of
the same story. It does the same thing (it does *more* things, even!), often in
ways that are a little different, because we have written thousands of line of
code using `EcologicalNetworks` and that gave us a really good understanding of
how *not* to write this package. So this is the second attempt.
