```@meta
DocTestSetup = quote
    using Serendip
end
```

# Introduction to Serendip

Serendip is a Finite Element library written in Julia language. The purpose of this library is to aid the research of new algorithms for the finite element method. Currently this library solves static and dynamic analyses in two and three-dimensions.

## Installation and basic usage

Install the package using the package manager, type `]` and then:

```
] add Serendip
```

To use Serendip, type in Julia REPL:

```
using Serendip
```

To test, type `]` and then:

```
] test Serendip
```

## Development version

To install the development version, type `]` and then

```
add Serendip#main
```