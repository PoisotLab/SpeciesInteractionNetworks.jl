# Random draws

!!! abstract

    The methods in this page draw random binary networks from probabilistic networks

```@docs
randomdraws
```

## Null models for hypothesis testing

The usual null models of network structure can be generated using a
[`PermutationConstraint`](@ref) as the first argument. Although technically
speaking, these are not permutations, it is nevertheless useful to map the
constraints from one method of network generation to another.

```@docs
nullmodel
```

## Linear filtering

The linear filter from [Stock2017Linear](@citet) is a general formulation of the
null models presented above.

```@docs
linearfilter
```
