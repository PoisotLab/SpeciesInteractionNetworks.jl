# Species Interaction Networks

The `SpeciesInteractionNetworks` package enables analyses of species interaction
networks in *Julia*. It is a re-master (director's cut?) of
[`EcologicalNetworks`](https://github.com/PoisotLab/EcologicalNetworks.jl/).

!!! question "Why a new package?"

    Code rot, experience, and truth in advertising. `EcologicalNetworks` was initially written in the days of *Julia* 0.5, and never really caught up with changes in the language (or with accumulated experience in writing code that is easy to maintain). After realizing the the former's type system was often preventing elegant dispatch and was due for a re-write, it made sense to initiate a susbtantial refactor. Additionally, there are ecological networks that do not describe species interactions, which were not covered by the package.

The package is built around a type system for species interaction networks,
which is intended to capture the different types of data and communities
ecologists need to handle. This makes the package extensible, both by writing
additional methods with a very fine-tuned dispatch, or by adding additional
types that should work out of the box (or be very close to).

This package is a library for the analysis of ecological networks. On purpose,
we do not provide "wrapper"-type functions that would perform an entire
analysis. We experimented with this idea during development, and rapidly
realized that even for the most simple research project, we needed to make small
tweaks that made the wrappers a nuisance. We decided to give you the parts, and
it's your job to build the kick-ass spaceship.

The measures in the documentation are organized by level of organisation,
following the convention of micro (species-level), meso (involving more than one
species), macro (network level) and meta (involving mutliple networks).