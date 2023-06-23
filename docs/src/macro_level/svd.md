# SVD and complexity

!!! abstract

    When generating networks, in particular when using [`randomdraws`](@ref) following [`nullmodel`](@ref), there is a chance that some species will have no interactions. In some cases, it may be relevant to identify and remove these species. The methods presented in this page offer a way to do this.

## Detection of disconnected species

```@docs
rdpg
complexity
```
