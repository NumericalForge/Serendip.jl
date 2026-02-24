# Getting Started

Serendip is a finite element library for research and engineering prototyping in Julia.

## Install

```julia
using Pkg
Pkg.add("Serendip")
```

For local development in this repository:

```julia
using Pkg
Pkg.develop(path=".")
```

## Load package

```julia
using Serendip
```

## First checks

```julia
Pkg.test("Serendip")
```

## Documentation scope

This documentation is currently organized in two tracks:

1. A curated manual with workflows and practical guidance.
2. A full API reference generated from exported symbols and docstrings.
