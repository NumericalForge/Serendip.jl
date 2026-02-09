<!-- # Mesh generation

## Node

### Node struct

```@docs
Node
```

### Node constructors

```@docs
Node()
Node(::Real,::Real,::Real)
Node(::AbstractArray)
```

### Node functions

```@docs
copy(::Node)
tag!(::Node, ::String)
tag!(::Vector{Node}, ::String)
add_dof
```

## Cell

### Cell struct

```@docs
Cell
```

### Cell constructors

```@docs
Cell(::CellShape, ::Vector{Node})
```

### Cell functions

```@docs
copy(::Cell)
tag!(::Cell, ::String)
tag!(::Vector{Cell}, ::String)
```

## Blocks

### Block struct
```@docs
Block
```

### Block constructors
```@docs
Block(::Array{Real})
```


### Block functions
```@docs
copy(::Block)
tag!(::Block, ::String)
tag!(::Vector{Block}, ::String)
array(::Block)
mirror(::Block)
polar(::Block)
rotate(::Block)
scale(::Block)
extrude(::Block)
```

## Mesh

### Mesh struct

```@docs
Mesh
```

### Mesh constructors

```@docs
Mesh(::Array{Real}, ::Array{Vector{Int64},1}, ::Vector{CellShape})
Mesh(::Block)
```

### Mesh functions -->