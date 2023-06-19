# Structural models

!!! abstract

    These functions generate random networks based on structural models. The models are given as types (specifically subtypes of [`StructuralModel`](@ref)), and return a network the type of which depends on the generative algorithm, when passed to the [`structuralmodel`](@ref) function.

## Models

```@docs
StructuralModel
NicheModel
```

## The generative function

```@docs
structuralmodel
```
