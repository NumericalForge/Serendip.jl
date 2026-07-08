# Serendip - Finite element code

<!-- <p align="center">
  <img src="docs/src/assets/Serendip2.png" />
</p> -->

Serendip is a Finite Element library written in Julia language. The purpose of this library is to aid the research of new algorithms for the finite element method. Currently, this library performs static and dynamic analyses in two and three dimensions.

Generic chart plotting in Serendip is provided by `QuickCharts` and re-exported from `Serendip`. Domain-specific visualization remains in Serendip through `DomainPlot`.

## Installation

Serendip currently requires Julia `1.12` or newer.

```julia
using Pkg
Pkg.add("Serendip")
```

<!-- [Documentation](https://NumericalForge.github.io/Serendip.jl/dev/) -->
