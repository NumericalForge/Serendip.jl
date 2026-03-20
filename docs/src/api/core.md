```@meta
CurrentModule = Serendip
```

# Core Structures

## Module Helpers

- `@t_str`

- `@L_str`

- `@run_files`

## Basic FEM Objects

```@raw html
<details>
<summary><code>Context</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Context
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>Constitutive</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Constitutive
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>Dof</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Dof
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>Node</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Node
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>Ip</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Ip
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>Element</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Element
```

```@raw html
</details>
```


## Node and Element Utilities

```@raw html
<details>
<summary><code>add_dof</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_dof
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>get_dof</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.get_dof
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>get_values</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.get_values
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>change_quadrature</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.change_quadrature
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>get_ips</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.get_ips
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>nearest</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.nearest
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>get_coords</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.get_coords
```

```@raw html
</details>
```


## Model Assembly

```@raw html
<details>
<summary><code>RegionMapper</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.RegionMapper
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>RegionModel</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.RegionModel
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>Stage</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.Stage
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>FEModel</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.FEModel
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>SolverSettings</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.SolverSettings
```

```@raw html
</details>
```


## Workflow Helpers

```@raw html
<details>
<summary><code>add_mapping</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_mapping
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>add_map</code></summary>
```

```@autodocs; canonical=false
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_map
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>add_stage</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_stage
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>add_bc</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_bc
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>add_logger</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_logger
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>add_monitor</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.add_monitor
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>run</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.run
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>reset_displacements</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.reset_displacements
```

```@raw html
</details>
```


```@raw html
<details>
<summary><code>save</code></summary>
```

```@autodocs
Modules = [Serendip]
Private = false
Public = true
Filter = x -> x === Serendip.save
```

```@raw html
</details>
```

