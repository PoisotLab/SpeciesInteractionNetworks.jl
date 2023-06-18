# Species Interaction Networks

The `SpeciesInteractionNetworks` package enables analyses of species interaction
networks in *Julia*. The list of implemented measures (and more broadly, the
design of the package) closely follows the recommendations in
[Delmas2019Analysing](@citet).

The measures in the documentation are organized by level of organisation,
following a simple convention. *Micro*-level measures return information at the
level of the node. Typically, these will return a value for each species in the
network. *Meso*-level measures return information about more than one species,
but not about the entire network. Examples are shortest paths, motif
enumeration, etc. *Macro*-level measures are summaries of the entire network,
and typically return one value for the entire network. Finally, *meta*-level
measures are intended to be applied to a collection of networks; network
dissimilarity measures are one example.

This package is a library for the analysis of ecological networks. On purpose,
we do not provide "wrapper"-type functions that would perform an entire
analysis. We experimented with this idea during development, and rapidly
realized that even for the most simple research project, we needed to make small
tweaks that made the wrappers a nuisance. We decided to give you the parts, and
it's your job to build the kick-ass spaceship.

The package is built around a type system for species interaction networks,
which is intended to capture the different types of data and communities
ecologists need to handle. This makes the package extensible, both by writing
additional methods with a very fine-tuned dispatch, or by adding additional
types that should work out of the box (or be very close to).

!!! question "Why a new package?"

    This package is a re-master (director's cut?) of [`EcologicalNetworks`](https://github.com/PoisotLab/EcologicalNetworks.jl/). Why? Code rot, experience, and truth in advertising. `EcologicalNetworks` was initially written in the days of *Julia* 0.5, and never really caught up with changes in the language (or with accumulated experience in writing code that is easy to maintain). After realizing the the former's type system was often preventing elegant dispatch and was due for a re-write, it made sense to initiate a susbtantial refactor. Additionally, there are ecological networks that do not describe species interactions, which were not covered by the package.
