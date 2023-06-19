# Measures of centrality

!!! abstract

    Centrality can help in quantifying the importance of species in a network. These function will measure the centrality of all species under different algorithms. There is a single wrapper function called `centrality`, which uses an optional first argument to specify the algorithm to use.

The centrality scores are returned so that they *sum* to one. This is intended
to make sure that *within a network*, the values for different nodes are
comparable.

## Implemented algorithms

```@docs
CentralityMethod
```

```@docs
KatzCentrality
EigenvectorCentrality
ClosenessCentrality
ResidualClosenessCentrality
GeneralizedClosenessCentrality
```

## Centrality function

```@docs
centrality
```
