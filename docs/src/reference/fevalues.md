```@meta
CurrentModule = Ferrite
DocTestSetup = :(using Ferrite)
```

# FEValues

## Main types
[`CellValues`](@ref) and [`FacetValues`](@ref) are the most common
subtypes of `Ferrite.AbstractValues`. For more details about how
these work, please see the related [topic guide](@ref fevalues_topicguide).

```@docs
CellValues
FacetValues
```

!!! warning "Embedded API"
    Currently, embedded `FEValues` returns `SArray`s, which behave differently
    from the `Tensor`s for normal value. In the future, we expect to return
    an `AbstractTensor`, this change may happen in a minor release, and the
    API for embedded `FEValues` should therefore be considered experimental.

## Applicable functions
The following functions are applicable to both `CellValues`
and `FacetValues`.

```@docs
reinit!
getnquadpoints
getdetJdV

shape_value(::Ferrite.AbstractValues, ::Int, ::Int)
shape_gradient(::Ferrite.AbstractValues, ::Int, ::Int)
shape_symmetric_gradient
shape_divergence
shape_curl
geometric_value

function_value
function_gradient
function_symmetric_gradient
function_divergence
function_curl
spatial_coordinate
```

In addition, there are some methods that are unique for `FacetValues`.

```@docs
Ferrite.getcurrentfacet
getnormal
```

## [InterfaceValues](@id reference-interfacevalues)

All of the methods for [`FacetValues`](@ref) apply for `InterfaceValues` as well.
In addition, there are some methods that are unique for `InterfaceValues`:

```@docs
InterfaceValues
shape_value_average
shape_value_jump
shape_gradient_average
shape_gradient_jump
function_value_average
function_value_jump
function_gradient_average
function_gradient_jump
```
