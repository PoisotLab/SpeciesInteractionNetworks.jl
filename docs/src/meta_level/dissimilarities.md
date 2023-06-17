# List of dissimilarity functions

!!! abstract

    The methods in this page are from [Koleff2003Measuring](@citet); when called on the output of [`betadiversity`](@ref), they will return a (dis)similarity score.

In the original paper [Koleff2003Measuring](@cite), the measures are expressed
as a function of ``a`` (the number of shared elements), ``b`` (the number of
elements unique to the right/neighbouring set), and ``c`` (the number of
elements unique to the left/focal set). In this documentation, we use a slightly
different convention, where ``c = L`` (left), ``a = S`` (shared), and ``b = R``
(right).

## List of measures

```@docs
KGL01
KGL08
```
